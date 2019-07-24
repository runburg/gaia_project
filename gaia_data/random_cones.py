#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generates random cones.

Author: Jack Runburg
Date: 21-06-2019 12:45


"""

import numpy as np
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import astropy.io.ascii as ascii
from astropy.table import Table, Column

customSimbad = Simbad()
customSimbad.add_votable_fields('morphtype', 'dim_majaxis', 'dim_minaxis', 'distance', 'pmdec', 'pmra')
# print(customSimbad.get_votable_fields())


def stringdms_to_deg(list):
    """Convert from dms to deg given a list of the three componenets."""
    dd = []
    for item in list:
        dd.append(float(item.split()[0]) + float(item.split()[1])/60 + float(item.split()[2])/3600)
    return dd


def stringhms_to_deg(list):
    """Convert from dms to deg given a list of the three componenets."""
    dd = []
    for item in list:
        dd.append((float(item.split()[0]) + float(item.split()[1])/60 + float(item.split()[2])/3600)*15)
    return dd


def pm_mag(pmdec, pmra):
    """Calculate simple magnitude estimate of proper motion; delete Nan rows as well."""
    ra_index = ~np.isnan(pmra)
    pmdec = pmdec[ra_index]
    pmra = pmra[ra_index]
    dec_index = ~np.isnan(pmdec)
    pmdec = pmdec[dec_index]
    pmra = pmra[dec_index]

    return np.array([np.sqrt(dec**2 + ra**2) for (dec, ra) in zip(pmdec, pmra)])


def reject_outliers(values, n=5):
    """Reject the outliers of an array that are n sigmas away from the mean."""
    values = values[~np.isnan(values)]
    return values[abs(values - np.mean(values)) < n * np.std(values)]


def gaia_search(simbad_table):
    """Given a table from Simbad, return gaia cone search around object."""
    r = []
    for object in simbad_table:
        coords = coord.SkyCoord(ra=object['RA'], dec=object['DEC'], unit=(u.degree, u.degree), frame='icrs')
        # radius = u.Quantity(object['GALDIM_MAJAXIS']/2, u.arcmin).to(u.degree)
        radius = 0.5 * u.degree
        j = Gaia.cone_search_async(coords, radius)
        r.append(j.get_results())
    return r


def check_if_in_galactic_plane(ra, dec):
    """Check if in plane and return new coordinates if not."""
    c_icrs = coord.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    if c_icrs.galactic.b < 2.5 and c_icrs.galactic.b > -2.5:
        dec += 10
    return ra, dec


# generate random cones
rando = np.random.rand(8, 3)
rando[:, 0] = np.interp(rando[:, 0], (rando[:, 0].min(), rando[:, 0].max()), (0, 360))
rando[:, 1] = np.interp(rando[:, 1], (rando[:, 1].min(), rando[:, 1].max()), (-90, 90))
rando[:, 2] = np.interp(rando[:, 2], (rando[:, 2].min(), rando[:, 2].max()), (4, 6))

random_table = Table()
random_table['RA'] = Column(rando[:, 0], unit='degree', description='random RA')
random_table['DEC'] = Column(rando[:, 1], unit='degree', description='random DEC')
random_table['GALDIM_MAJAXIS'] = Column(rando[:, 2], unit='arcmin', description='random radius')

random_table['GAIA'] = gaia_search(random_table)
print("done querying ii")
random_table['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in random_table['GAIA']])

ascii.write(random_table, output='./randomcone_info.ecsv', format='ecsv', overwrite=True)

with open('randomcone_vels.npz', 'wb') as outfile:
    np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in random_table['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in random_table['GAIA']])
