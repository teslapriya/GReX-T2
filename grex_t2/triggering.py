import os.path
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import logging as logger

logger.basicConfig(filename="logs/output.log", encoding="utf-8", level=logger.DEBUG)

# TODO should move to calib constants
Ddish = 4.5
dmax = 273.0
freq = 1428.0
lamb = 3e2 / freq
nbeam = 256
beam_separation = 1  # arcminutes
HA_pointing = 0.0
npx = 1000
theta_min = -(nbeam / 2.0) * beam_separation


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

    if catalog is not None:
        assert os.path.exists(catalog), f"catalog file ({catalog}) not found"

        with open(catalog, "r") as reader:
            lines = reader.readlines()

        for i in np.arange(1, len(lines)):
            try:
                ra, dec, minsnr = lines[i].split()
                c = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame="icrs")
                coords.append(c)
                snrs.append(float(minsnr))
            except:
                print("")
    else:
        logger.warning("No catalog found. Will not filter output based on a catalog.")
        print("No catalog found. Will not filter output based on a catalog.")

    return coords, snrs


def check_clustered_sources(tab, coords, snrs, beam_model=None, do_check=True):
    """Reduces tab according to sources

    Parameters
    ----------
    tab: astropy table of clustered candidates
    coords: list of astropy coordinates of sources
    snrs: list of snr flags for sources (-1 is inf, -2 is primary beam flagging)
    beam_model: standard beam model input to beams_coord (default None)

    Returns
    -------
    filtered astropy table of candidates
    """

    if do_check is False:
        return tab

    good = [True] * len(tab)

    # mjds and ibeam and snr are the three columns of interest
    mjd = tab["mjds"]
    ibeam = tab["ibeam"]
    snr = tab["snr"]
    ncand = len(mjd)
    is_not_src = [True] * ncand

    # select based on beam and snr
    for i in np.arange(ncand):
        for j in np.arange(len(snrs)):
            if snrs[j] == -2.0:
                beams, resps = primary_beams_coord(
                    coords[j].ra.deg,
                    coords[j].dec.deg,
                    mjd[i],
                    beam_model=beam_model,
                )
            else:
                beams, resps = beams_coord(
                    coords[j].ra.deg,
                    coords[j].dec.deg,
                    mjd[i],
                    beam_model=beam_model,
                )
            if ibeam[i] + 1 in beams:  # +1 is for model v T1/T2 offset
                if snr[i] < snrs[j]:
                    is_not_src[i] = False
                if snrs[j] == -1.0:
                    is_not_src[i] = False
                if snrs[j] == -2.0:
                    is_not_src[i] = False

    tab_out = tab[is_not_src]

    print(f"Filtering from {len(tab)} to {len(tab_out)} candidates (source check).")
    return tab_out
