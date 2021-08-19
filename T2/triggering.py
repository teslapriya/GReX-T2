#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import json
import os.path
import numpy as np
#from sklearn import cluster  # for dbscan
import hdbscan
from astropy import time, coordinates
from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError
from astropy.coordinates import SkyCoord, ITRS, EarthLocation
from astropy.time import Time
from astropy import units as u
from astropy import wcs
from dsautils import dsa_store
from progress.bar import Bar
import dsacalib.constants as ct
import dsautils.cnf as cnf
from T2 import cluster_heimdall
ds = dsa_store.DsaStore()
import dsautils.dsa_syslog as dsl
logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')
MY_CNF = cnf.Conf()
CORR_CNF = MY_CNF.get('corr')

# TODO should move to calib constants
Ddish = 4.5
dmax=273.
freq=1428.
lamb = 3e2 / freq
nbeam = 256
beam_separation = 1 # arcminutes	
HA_pointing = 0.0
npx = 1000
theta_min = -(nbeam/2.)*beam_separation

def parse_catalog(catalog):
    """Parse source catalog

    Parameters
    ----------
    catalog: string path to file
        Needs to have format <ra_sexagesimal> <dec_sexagesimal> <SNR_flag>

    Returns
    -------
    list of astropy coordinates, list of SNRs
    """
    
    coords = []
    snrs = []

    with open(catalog,"r") as reader:
        lines = reader.readlines()

    for i in np.arange(1,len(lines)):

        try:
            ra,dec,minsnr = lines[i].split()
            c = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
            coords.append(c)
            snrs.append(float(minsnr))
        except:
            print("")

    return coords,snrs
    

def get_pointing_declination(tol=0.25):
    """Gets the pointing declination from the commanded antenna elevations.

    Parameters
    ----------
    tol : float
        The tolerance for discrepancies in the antenna pointing and commanded
        elevations, in degrees.

    Returns
    -------
    astropy quantity
        The pointing declination, in degrees or equivalent.
    """
    commanded_els = np.zeros(len(CORR_CNF['antenna_order']))
    for idx, ant in CORR_CNF['antenna_order'].items():
        try:
            antmc = ds.get_dict('/mon/ant/{0}'.format(ant))
            a1 = np.abs(antmc['ant_el'] - antmc['ant_cmd_el'])
        except:
            a1 = 2.*tol
        if a1 < tol:
            commanded_els[idx] = antmc['ant_cmd_el']
        else:
            commanded_els[idx] = np.nan

    pt_el = np.nanmedian(commanded_els)
    if pt_el is not np.nan:
        pt_dec = ct.OVRO_LAT*u.rad + pt_el*u.deg - 90*u.deg
    else:
        pt_el = CORR_CNF['pt_dec']
    return pt_dec


# create WCS
def create_WCS(c,cdelt_deg):
    """Creates astropy wcs object for beams
    
    Parameters
    ----------
    c : astropy SkyCoord object
        Center coordinates
    cdelt_deg: float
        Separation between beams in degrees

    Returns
    -------
    astropy wcs object
        wcs corresponding to image of beams.
    """

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [npx/2, npx/2]
    w.wcs.cdelt = np.array([cdelt_deg, cdelt_deg])
    w.wcs.crval = [c.ra.deg, c.dec.deg]
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]

    return w

def gauss1d(x,mu,sigma):
    """1D Gaussian envelope
    
    Parameters
    ----------
    x: float
        x-value to evaluate at
    mu: float
        mean
    sigma: float
        sigma

    Returns
    -------
    float
        value of 1d gaussian
    """
    
    return np.exp(-(x-mu)**2 / (2.*sigma**2))

def gaussian2D(coords,  # x and y coordinates for each image.
               amplitude=1,  # Highest intensity in image.
               xo=0,  # x-coordinate of peak centre.
               yo=0,  # y-coordinate of peak centre.
               sigma_x=1,  # Standard deviation in x.
               sigma_y=1,  # Standard deviation in y.
               rho=0,  # Correlation coefficient.
               offset=0,
               rot=0):  # rotation in degrees.

    """Returns 2D Gaussian envelope 
    
    Parameters
    ----------
    coords : (x,y) tuple of meshgrid coordinates for each image
    amplitude: float
        Highest intensity in image. (default 1.0)
    x0: int
        x-coordinate of peak centre. (default 0)
    y0: int
        y-coordinate of peak centre.  (default 0)
    sigma_x: float
        Standard deviation in x. (default 1.0)
    sigma_y: float
        Standard deviation in y. (default 1.0)
    rho: float
        Correlation coefficient between x and y, stretches Gaussian (default 0)
    offset: float
        Baseline offset (default 0.0)
    rot: float
        Rotation of Gaussian in degrees (default 0.0)

    Returns
    -------
    Squeezed 2D grid representing Gaussian.

    """
    
    x, y = coords

    rot = np.deg2rad(rot)

    x_ = np.cos(rot)*x - y*np.sin(rot)
    y_ = np.sin(rot)*x + np.cos(rot)*y

    xo = float(xo)
    yo = float(yo)

    xo_ = np.cos(rot)*xo - yo*np.sin(rot) 
    yo_ = np.sin(rot)*xo + np.cos(rot)*yo

    x,y,xo,yo = x_,y_,xo_,yo_

    # Create covariance matrix
    mat_cov = [[sigma_x**2, rho * sigma_x * sigma_y],
               [rho * sigma_x * sigma_y, sigma_y**2]]
    mat_cov = np.asarray(mat_cov)
    # Find its inverse
    mat_cov_inv = np.linalg.inv(mat_cov)

    # PB We stack the coordinates along the last axis
    mat_coords = np.stack((x - xo, y - yo), axis=-1)

    G = amplitude * np.exp(-0.5*np.matmul(np.matmul(mat_coords[:, :, np.newaxis, :],mat_cov_inv),mat_coords[..., np.newaxis])) + offset
    return G.squeeze()


# this is in MD/alt coordinates
def get_2Dbeam_model():
    """Returns a 2D beam model image composed of 2D Gaussians

    Parameters
    ----------

    Returns
    -------
    numpy array
        [nbeam, npix, npix] array of beam response values
    """


    # Arcminutes 
    sb_width_fwhm = lamb / dmax * 180 / np.pi * 120

    # Arcminutes
    primary_width_fwhm = lamb / Ddish * 180 / np.pi * 60

    beam_formed_min, beam_formed_max = -2.*primary_width_fwhm, 2.*primary_width_fwhm

    theta = np.linspace(beam_formed_min, beam_formed_max, npx)

    coords = np.meshgrid(theta,theta)

    G = gaussian2D(coords, xo=0, yo=0, sigma_x=sb_width_fwhm/2.355, sigma_y=primary_width_fwhm/2.355)
    beam_env = gauss1d(theta, 0, primary_width_fwhm/2.355)
    mus = (np.arange(256)*beam_separation + theta_min)
    beam_val = np.zeros([nbeam, len(theta), len(theta)])
    Genv = gaussian2D(coords, xo=0, yo=0, sigma_x=primary_width_fwhm/2.355, sigma_y=primary_width_fwhm/2.355)

    with Bar('Calculating beam model...',suffix='%(percent).1f%% - %(eta)ds',max=len(mus)) as bar:
        for ii, mu in enumerate(mus):
            G = gaussian2D(coords, xo=mu, yo=0, sigma_x=sb_width_fwhm/2.355, sigma_y=primary_width_fwhm/2.355*10)
            #beam_val[ii] = G*Genv
            beam_val[ii] = G
            bar.next()

    return beam_val 

def beams_coord(ra_deg,dec_deg,mjd,dec=None,response=0.1,beam_model=None):
    """Provides a list of beams with response above a given value at a coordinate

    Parameters
    ----------
    ra_deg: float (degrees)
        RA coordinate (ICRS) in degrees
    dec_deg: float (degrees)
        DEC coordinate (ICRS) in degrees
    mjd: float
        mjd to calculate beams at 
    dec: float (degrees)
        Pointing declination (default is to try to figure it out from etcd/cnf)
    response: float
        Response limit above which beam is returned, as fraction of peak (default 0.1)
    beam_model: output of get_2Dbeam_model (default None)
    
        
    Returns
    -------
    list
        list of beams that contribute at position (can be empty)
    """

    c = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    RA_radians = np.deg2rad(c.ra.value)    

    if dec is None:
        Dec = get_pointing_declination()
        Dec = Dec.value*180./np.pi
    else:
        Dec = dec

    
    # Arcminutes
    sb_width_fwhm = lamb / dmax * 180 / np.pi * 120

    # Arcminutes
    primary_width_fwhm = lamb / Ddish * 180 / np.pi * 60
    
    beam_formed_min, beam_formed_max = -2.*primary_width_fwhm, 2.*primary_width_fwhm
    theta = np.linspace(beam_formed_min, beam_formed_max, npx)
    coords = np.meshgrid(theta,theta)
        
    dalt_arr, dMD_arr = coords[0], coords[1] # these are alt, MD
    if beam_model is None:
        beam_val = get_2Dbeam_model() # this is in alt / MD, not HA
    else:
        try:
            beam_val = beam_model
        except:
            beam_val = get_2Dbeam_model() # this is in alt / MD, not HA
            
    # get results
    t = Time(mjd, format='mjd', scale='utc')        
    c_ITRS = c.transform_to(ITRS(obstime=t))
    local_ha = ct.OVRO_LON*u.rad - c_ITRS.spherical.lon
    HA_src = local_ha.deg+360.
    RA_pt =  (t.sidereal_time('apparent', longitude=ct.OVRO_LON*(180./np.pi)*u.deg)).deg
    coord_pt = SkyCoord(ra=RA_pt*u.deg, dec=Dec*u.deg, frame='icrs')
    w = create_WCS(coord_pt,(theta[1]-theta[0])/60.)
    xy = w.wcs_world2pix((RA_pt+HA_src), c.dec.deg, 0)

    if xy[1]>=npx:
        return [],[]
    if xy[1]<0:
        return [],[]
    if xy[0]>=npx:
        return [],[]
    if xy[0]<0:
        return [],[]
    
    result = beam_val[:, int(xy[1]), int(xy[0])]
    beams_out = []
    result = result.ravel()
    vals = []
    i = 0
    for val in result:
        if val>response:
            beams_out.append(i)
            vals.append(val)
        i += 1

    return beams_out,vals


def check_clustered_sources(tab,coords,snrs,beam_model=None):
    """ Reduces tab according to sources

    Parameters
    ----------
    tab: astropy table of clustered candidates
    coords: list of astropy coordinates of sources
    snrs: list of snr flags for sources (-1 is inf)
    beam_model: standard beam model input to beams_coord (default None)

    Returns
    -------
    filtered astropy table of candidates
    """

    good = [True] * len(tab)
    
    # mjds and ibeam and snr are the three columns of interest
    mjd = tab['mjds']
    ibeam = tab['ibeam']
    snr = tab['snr']
    ncand = len(mjd)
    is_not_src = [True] * ncand

    # select based on beam and snr
    for i in np.arange(ncand):

        for j in np.arange(len(snrs)):        
            beams,resps = beams_coord(coords[j].ra.deg,coords[j].dec.deg,mjd[i],beam_model=beam_model)
            if ibeam[i]+1 in beams: # +1 is for model v T1/T2 offset
                if snr[i]<snrs[j]:
                    is_not_src[i] = False
                if snrs[j]==-1.:
                    is_not_src[i] = False

    tab_out = tab[is_not_src]

    print(f'Filtering from {len(tab)} to {len(tab_out)} candidates (source check).')
    return tab_out
                    
def test_on_clustered(T2_cands_file,source_check=False,catalog=None,beam_model=None,max_ncl=10,outputfile='tmp.dat'):

    tab = cluster_heimdall.parse_candsfile(T2_cands_file)
    tab2 = cluster_heimdall.filter_clustered(tab,min_snr=8.0,min_dm=20.0,max_ibox=33.0)
    tab2['gulp_ncl'] = [len(tab)] * len(tab2)

    if len(tab2) and len(tab2)<max_ncl:
        
        maxsnr = tab2['snr'].max()
        imaxsnr = np.where(tab2['snr'] == maxsnr)[0][0]
        tab3 = tab2[imaxsnr:imaxsnr+1]

        if source_check is True and catalog is not None:

            coords,snrs = parse_catalog(catalog)
            tab4 = check_clustered_sources(tab3,coords,snrs,beam_model=beam_model)

            if len(tab4):
                tab4.write(outputfile, format='ascii.no_header')

        else:
            tab3.write(outputfile, format='ascii.no_header')

            
