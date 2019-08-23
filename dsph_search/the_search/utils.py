#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for the search.

Author: Jack Runburg
Date: 22-08-2019 14:30


"""

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import numpy as np
from random import random as random

def gaia_search(ra, dec, radius=0.5):
    """Given a table from Simbad, return gaia cone search around object."""
    # radius = radius * u.degree
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    # radius = u.Quantity(object['GALDIM_MAJAXIS']/2, u.arcmin).to(u.degree)
    job = Gaia.launch_job_async(f"SELECT TOP 50000 \
                                gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec, \
                                gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error, \
                                gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error, \
                                gaia_source.bp_rp \
                                FROM gaiadr2.gaia_source \
                                WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{ra},{dec},{radius}))=1")

    return job


def random_cones_outside_galactic_plane(limit=15):
    """Check if in plane and return new coordinates if not."""
    # galactic longitude
    l = random() * 360
    # galactic latitude
    b = (random() - 0.5) * 90

    # check that no galactic latitudes are within limit deg of galactic plane
    if abs(b) < limit:
        b += limit * np.sign(b)

    c_gal = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
    icrs_coords = (c_gal.icrs.ra.value, c_gal.icrs.dec.value)

    return icrs_coords


if __name__ == '__main__':
    ra, dec = random_cones_outside_galactic_plane()
    print(ra, dec)
