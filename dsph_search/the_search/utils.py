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
import warnings


def gaia_search(ra, dec, name, radius=0.5, sigma=5, pm_threshold=5, dump_to_file=True):
    """Given a table from Simbad, return gaia cone search around object."""
    warnings.filterwarnings("ignore", module='astropy.*')
    # radius = radius * u.degree
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    # radius = u.Quantity(object['GALDIM_MAJAXIS']/2, u.arcmin).to(u.degree)
    job = Gaia.launch_job_async(f"SELECT TOP 50000 \
                                gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec, \
                                gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error, \
                                gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error, \
                                gaia_source.bp_rp \
                                FROM gaiadr2.gaia_source \
                                WHERE \
                                CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{ra},{dec},{radius}))=1 AND  (gaiadr2.gaia_source.parallax - gaiadr2.gaia_source.parallax_error * {sigma} <= 0) AND (SQRT(POWER(gaiadr2.gaia_source.pmra, 2) + POWER(gaiadr2.gaia_source.pmdec, 2)) <= {pm_threshold})", dump_to_file=dump_to_file, output_file=f'./candidates/{name}/vots/{name}_{round(radius*100)}.vot')

    return job


def random_cones_outside_galactic_plane(limit=15):
    """Check if in plane and return new coordinates if not."""
    # galactic longitude
    l = random() * 360
    # galactic latitude
    b = (random() - 0.5) * (90-limit)
    # check that no galactic latitudes are within limit deg of galactic plane

    if b<0:
        b -= limit
    else:
        b += limit

    c_gal = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
    icrs_coords = (c_gal.icrs.ra.value, c_gal.icrs.dec.value)

    return icrs_coords


def unmask(data):
    """Return an unmasked table for cuts."""
    data = data[[~obj.mask for obj in data]]

    return data


if __name__ == '__main__':
    ra, dec = random_cones_outside_galactic_plane()
    gaia_search(90, 90, 'test', dump_to_file=False)
    print(ra, dec)
