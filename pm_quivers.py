#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make histograms of proper motion.

Author: Jack Runburg
Date: 02-08-2019 11:14


"""

import gaia_functions as g
import astropy.io.ascii as ascii
import numpy as np
import matplotlib.pyplot as plt

def main():
    dwarf_specs = ascii.read('./gaia_data/dwarf_info.ecsv', format='ecsv')
    xcoords = dwarf_specs['RA']
    ycoords = dwarf_specs['DEC']
    titles = dwarf_specs['MAIN_ID']

    # load pm values
    with np.load('./gaia_data/dwarf_vels.npz', allow_pickle=True) as infile:
        # vels = infile['vels']
        pmra = infile['pmra']
        pmdec = infile['pmdec']
        ra = infile['ra']
        dec = infile['dec']

    cuts = [0.1, 0.5, 1, 5, 10]
    for cut in cuts:
        cut_name = dict(zip(cuts, ['tenth', 'half', '1', '5', '10']))
        g.pm_vector_plot([ra[0]], [dec[0]], [pmra[0]], [pmdec[0]], f'./plots/draco_quiver_{cut_name[cut]}.pdf', names=titles, cut=cut)



if __name__ == "__main__":
    main()
