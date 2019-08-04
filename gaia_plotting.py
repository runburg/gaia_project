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
for radius in radii:
    # g.make_pm_maps(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/known_dwarf_histos_{radius}.pdf', num_cones=8, titles=1)
    # g.make_pm_maps(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/random_cone_histos_{radius}.pdf', num_cones=8)
    g.parallax_histograms(f'./gaia_data/dwarf_info_{radius}.ecsv', f'./gaia_data/dwarf_vels_{radius}.npz', f'./plots/known_dwarf_parallax_histos_{radius}.pdf', titles=1)
    g.parallax_histograms(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/random_cone_parallax_histos_{radius}.pdf')
