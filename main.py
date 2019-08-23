#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Query SIMBAD and GAIA for astrophysical data.

Author: Jack Runburg
Date: 22-06-2019 14:29


"""
import gaia_functions as g
import numpy as np
import astropy.io.ascii as ascii
from astroquery.simbad import Simbad

# ## known dwarfs
# dwarfs to look at
dwarf_list = ["Dra dSph", "Leo I dSph", "Leo II dSph", "UMi dSph", "Sculptor dSph", "Carina dSph", "Sextans dSph", "Fornax dSph", "Boo dSph", "Cetus II", "Col I", "Gru II", "Coma Dwarf Galaxy", "Segue I"]
# dwarf_list = ['Dra dSph', 'Leo I dSph']
# simbad data requested
customSimbad = Simbad()
customSimbad.add_votable_fields('morphtype', 'dim_majaxis', 'dim_minaxis', 'distance', 'pmdec', 'pmra')

# query simbad
dwarf_specs = customSimbad.query_objects(dwarf_list)

# change positions to degrees
print(dwarf_specs['MAIN_ID','RA', 'DEC'])

dwarf_specs['RA'] = g.stringhms_to_deg(dwarf_specs['RA'])
dwarf_specs['DEC'] = g.stringdms_to_deg(dwarf_specs['DEC'])

print(dwarf_specs['MAIN_ID','RA', 'DEC'])

radii = [0.5, 0.1, 0.05]
for radius in radii:
    # now query gaia
    dwarf_specs['GAIA'] = g.gaia_search(dwarf_specs, radius=radius)
    print("done querying i")


    # only include the objects gaia measured pms for
    dwarf_specs['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in dwarf_specs['GAIA']])
    # & (obj['parallax'] < 1/20)
    # print(dwarf_specs[np.less(dwarf_specs['GAIA'][0]['parallax'], 0.5)])

    # save pms and Simbad tables
    ascii.write(dwarf_specs, output=f'./gaia_data/dwarf_info_{radius}.ecsv', format='ecsv', overwrite=True)
    with open(f'./gaia_data/dwarf_vels_{radius}.npz', 'wb') as outfile:
        np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in dwarf_specs['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in dwarf_specs['GAIA']], ra=[np.array(obj['ra']) for obj in dwarf_specs['GAIA']], dec=[np.array(obj['dec']) for obj in dwarf_specs['GAIA']], parallax=[np.array(obj['parallax']) for obj in dwarf_specs['GAIA']], parallax_error=[np.array(obj['parallax_error']) for obj in dwarf_specs['GAIA']])

# ## random cones
# generate random cone coordinates
# dwarf_specs['RA'], dwarf_specs['DEC'] = g.random_cones_outside_galactic_plane(len(dwarf_list))
dwarf_specs['RA'] = np.array([ 65.47629119236987, 326.48154125527896, 72.7146753092771, 139.1067038659881, 346.25127522798704, 289.1081609393706, 257.0704958874846, 332.56218687772304, 247.96320615330532, 248.8958744426278, 295.0776709675192, 24.216171811340296, 328.3579688732271, 257.77631791021463])
dwarf_specs['DEC'] = np.array([ 25.257528676684768, 10.474170121853948, -74.84684844097445, -14.851805273629752, 34.65722523571136, -36.555318085426165, 12.406413156130805, -1.2835395220206993, -4.7979387679353325, -20.173861396876788, 48.189809814586845, 34.87627566749339, 27.280421885116603, 63.91444684764842])

#
# for radius in radii:
#     # query gaia
#     dwarf_specs['GAIA'] = g.gaia_search(dwarf_specs, radius)
#     print("done querying ii")
#
#     # only include objects gaia meausred pms for
#     dwarf_specs['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in dwarf_specs['GAIA']])
#
#     # save pms and random coordinates
#     ascii.write(dwarf_specs, output=f'./gaia_data/randomcone_info_{radius}.ecsv', format='ecsv', overwrite=True)
#     with open(f'./gaia_data/randomcone_vels_{radius}.npz', 'wb') as outfile:
#         np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in dwarf_specs['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in dwarf_specs['GAIA']], ra=[np.array(obj['ra']) for obj in dwarf_specs['GAIA']], dec=[np.array(obj['dec']) for obj in dwarf_specs['GAIA']], parallax=[np.array(obj['parallax']) for obj in dwarf_specs['GAIA']], parallax_error=[np.array(obj['parallax_error']) for obj in dwarf_specs['GAIA']])
#
# print(dwarf_specs['RA'])
# print(dwarf_specs['DEC'])
