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
# now query gaia
dwarf_specs['GAIA'] = g.gaia_search(dwarf_specs)
print("done querying i")
print(dwarf_specs['GAIA'][1])

# only include the objects gaia measured pms for
dwarf_specs['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in dwarf_specs['GAIA']])
# & (obj['parallax'] < 1/20)
# print(dwarf_specs[np.less(dwarf_specs['GAIA'][0]['parallax'], 0.5)])

# save pms and Simbad tables
ascii.write(dwarf_specs, output='./dwarf_info.ecsv', format='ecsv', overwrite=True)
with open('dwarf_vels.npz', 'wb') as outfile:
    np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in dwarf_specs['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in dwarf_specs['GAIA']], ra=[np.array(obj['ra']) for obj in dwarf_specs['GAIA']], dec=[np.array(obj['dec']) for obj in dwarf_specs['GAIA']])

# # ## random cones
# # generate random cone coordinates
# dwarf_specs['RA'], dwarf_specs['DEC'] = g.random_cones_outside_galactic_plane(len(dwarf_list))
#
# # query gaia
# dwarf_specs['GAIA'] = g.gaia_search(dwarf_specs)
# print("done querying ii")
#
# # only include objects gaia meausred pms for
# dwarf_specs['GAIA'] = np.array([obj[~obj['pmra'].mask & ~obj['pmdec'].mask] for obj in dwarf_specs['GAIA']])
#
# # save pms and random coordinates
# ascii.write(dwarf_specs, output='./randomcone_info.ecsv', format='ecsv', overwrite=True)
# with open('randomcone_vels.npz', 'wb') as outfile:
#     np.savez(outfile, pmra=[np.array(obj['pmra']) for obj in dwarf_specs['GAIA']], pmdec=[np.array(obj['pmdec']) for obj in dwarf_specs['GAIA']], ra=[np.array(obj['ra']) for obj in dwarf_specs['GAIA']], dec=[np.array(obj['dec']) for obj in dwarf_specs['GAIA']])
