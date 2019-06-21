#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Add up all the measured radial velocities for a given dwarf.

Author: Jack Runburg
Date: 29-05-2019 09:30

Half-light radii from https://doi.org/10.1111/j.1365-2966.2008.13739.x
"""

import numpy as np
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
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


dwarf_list = ["Dra dSph", "Leo I dSph", "Leo II dSph", "UMi dSph", "Sculptor dSph", "Carina dSph", "Sextans dSph", "Fornax dSph", "Boo dSph"]
dwarf_specs = customSimbad.query_objects(dwarf_list)
dwarf_specs['RA'] = stringhms_to_deg(dwarf_specs['RA'])
dwarf_specs['DEC'] = stringdms_to_deg(dwarf_specs['DEC'])
dwarf_specs = dwarf_specs[~dwarf_specs['PMRA'].mask & ~dwarf_specs['PMDEC'].mask]
dwarf_specs['PM_MAG'] = pm_mag(dwarf_specs['PMDEC'], dwarf_specs['PMRA'])
dwarf_specs['Distance_distance'] = dwarf_specs['Distance_distance'] * 1000 * u.kpc

dwarf_specs['GAIA'] = gaia_search(dwarf_specs)
dwarf_specs['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in dwarf_specs['GAIA']])
dwarf_specs['VELS'] = [np.array(reject_outliers(pm_mag(obj['pmdec'], obj['pmra']), n=4)) for obj in dwarf_specs['GAIA']]

ascii.write(dwarf_specs, output='./dwarf_info.ecsv', format='ecsv', overwrite=True)

with open('dwarf_vels.npz', 'wb') as outfile:
    np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in dwarf_specs['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in dwarf_specs['GAIA']])

# generate random cones
rando = np.random.rand(8, 3)
rando[:, 0] = np.interp(rando[:, 0], (rando[:, 0].min(), rando[:, 0].max()), (dwarf_specs['RA'].min(), dwarf_specs['RA'].max()))
rando[:, 1] = np.interp(rando[:, 1], (rando[:, 1].min(), rando[:, 1].max()), (dwarf_specs['DEC'].min(), dwarf_specs['DEC'].max()))
rando[:, 2] = np.interp(rando[:, 2], (rando[:, 2].min(), rando[:, 2].max()), (dwarf_specs['GALDIM_MAJAXIS'].min(), dwarf_specs['GALDIM_MAJAXIS'].max()))

random_table = Table()
random_table['RA'] = Column(rando[:, 0], unit='degree', description='random RA')
random_table['DEC'] = Column(rando[:, 1], unit='degree', description='random DEC')
random_table['GALDIM_MAJAXIS'] = Column(rando[:, 2], unit='arcmin', description='random radius')

random_table['GAIA'] = gaia_search(random_table)
random_table['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in random_table['GAIA']])

ascii.write(random_table, output='./randomcone_info.ecsv', format='ecsv', overwrite=True)

with open('randomcone_vels.npz', 'wb') as outfile:
    np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in random_table['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in random_table['GAIA']])
