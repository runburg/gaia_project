#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plotting script for Gaia project temperature maps.

Author: Jack Runburg
Date: 22-06-2019 14:41


"""
import gaia_functions as g
import matplotlib.pyplot as plt

g.make_pm_maps('./gaia_data/dwarf_info.ecsv', './gaia_data/dwarf_vels.npz', './plots/known_dwarf_histos.pdf', num_cones=8, titles=1)
g.make_pm_maps('./gaia_data/randomcone_info.ecsv', './gaia_data/randomcone_vels.npz', './plots/random_cone_histos.pdf', num_cones=8)
plt.show()
