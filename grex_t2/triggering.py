import os.path
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import logging as logger

logger.basicConfig(filename="output.log", encoding="utf-8", level=logger.DEBUG)

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
