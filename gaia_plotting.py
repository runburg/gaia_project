#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plotting script for Gaia project temperature maps.

Author: Jack Runburg
Date: 22-06-2019 14:41


"""
import gaia_functions as g
import matplotlib.pyplot as plt
import glob

radii = [0.5, 0.1, 0.05]
dwarf = 7
dwarf_dict = {
    0: 'draco',
    1: 'leoI',
    2: 'leoB',
    3: 'umi',
    4: 'sculptor',
    5: 'carina',
    6: 'sextans',
    7: 'fornax',
    8: 'bootes',
    9: 'cetus',
    10: 'col',
    11: 'gru',
    12: 'coma',
    13: 'segue'
    }
cuts=[1, 3, 5]
# radii = [0.1]
for radius in radii:
    g.make_pm_maps(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/pm_maps/dwarf_histos_{radius}.pdf', num_cones=8, titles=1)
    print('dwarf pm maps done')
    g.make_pm_maps(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/pm_maps/random_histos_{radius}.pdf', num_cones=8)
    print('random pm maps done')
    
    g.comparison_plot(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/comparisons/dwarf_parallax_histos_{radius}.pdf', titles=1)
    print('dwarf comparison done')
    g.comparison_plot(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/comparisons/random_parallax_histos_{radius}.pdf')
    print('random comparison done')

    for dwarf in dwarf_dict.keys():
        g.draco_comparison_plot(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/cuts/{dwarf_dict[dwarf]}_comparison_plot_{radius}.pdf', cuts=cuts, dwarf=dwarf, titles=1, radius=radius)
        print('dwarf cuts done', dwarf)
        g.draco_comparison_plot(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/cuts/random{dwarf+1}_comparison_plot_{radius}.pdf', cuts=cuts, dwarf=dwarf, radius=radius)
        print('random cuts done', dwarf)
