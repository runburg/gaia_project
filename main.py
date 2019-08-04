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
dwarf_specs['RA'] = g.stringhms_to_deg(dwarf_specs['RA'])
dwarf_specs['DEC'] = g.stringdms_to_deg(dwarf_specs['DEC'])

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
        np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in dwarf_specs['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in dwarf_specs['GAIA']], ra=[np.array(obj['ra']) for obj in dwarf_specs['GAIA']], dec=[np.array(obj['dec']) for obj in dwarf_specs['GAIA']], parallax=[np.array(obj['parallax']) for obj in dwarf_specs['GAIA']])

# ## random cones
# generate random cone coordinates
# dwarf_specs['RA'], dwarf_specs['DEC'] = g.random_cones_outside_galactic_plane(len(dwarf_list))
dwarf_specs['RA'] = np.array([84.57279212729698, 356.9910306983101, 188.3915695473758, 88.64349279201299, 45.62649691781021, 91.96720336392542, 0.8329131387295983, 186.1409551779438, 349.6709546202338, 278.01130769960116, 353.1985200164832, 20.913809787077344, 353.2045291740327, 199.71224834689642])
dwarf_specs['DEC'] = np.array([-59.94783760354259, 32.75788075774204, 16.72734338992764, 21.344532629253383, -14.961102434666808, 7.655718227382378, 40.20721677646744, 24.047301344703605, 13.671687100623563, -44.30547942050413, 17.140936852386147, 18.623002469314425, -20.900113392451228, 30.376185740855554])


for radius in radii:
    # query gaia
    dwarf_specs['GAIA'] = g.gaia_search(dwarf_specs, radius)
    print("done querying ii")

    # only include objects gaia meausred pms for
    dwarf_specs['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in dwarf_specs['GAIA']])

    # save pms and random coordinates
    ascii.write(dwarf_specs, output=f'./gaia_data/randomcone_info_{radius}.ecsv', format='ecsv', overwrite=True)
    with open(f'./gaia_data/randomcone_vels_{radius}.npz', 'wb') as outfile:
        np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in dwarf_specs['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in dwarf_specs['GAIA']], ra=[np.array(obj['ra']) for obj in dwarf_specs['GAIA']], dec=[np.array(obj['dec']) for obj in dwarf_specs['GAIA']], parallax=[np.array(obj['parallax']) for obj in dwarf_specs['GAIA']])
