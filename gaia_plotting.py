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
dwarf = 0
dwarf_dict = {
    0: 'draco',
    1: 'leoI',
    2: 'leoB',
    3: 'umi',
    4: 'sculptor'
    }
# radii = [0.1]
for radius in radii:
    g.make_pm_maps(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/known_dwarf_histos_{radius}.pdf', num_cones=8, titles=1)
    g.make_pm_maps(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/random_cone_histos_{radius}.pdf', num_cones=8)
    g.comparison_plot(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/known_dwarf_parallax_histos_{radius}.pdf', titles=1)
    g.comparison_plot(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/random_cone_parallax_histos_{radius}.pdf')
    g.draco_comparison_plot(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/{dwarf_dict[dwarf]}_comparison_plot_{radius}.pdf', cuts=[0.05, 0.1, 0.5, 1], dwarf=dwarf, titles=1, radius=radius)
    g.draco_comparison_plot(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/random{dwarf+1}_comparison_plot_{radius}.pdf', cuts=[0.05, 0.1, 0.5, 1], dwarf=dwarf, radius=radius)
