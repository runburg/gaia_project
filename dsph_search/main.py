#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main script for running the dSph search and controlling output.

Author: Jack Runburg
Date: 22-08-2019 14:03

"""

import numpy as np
import os, glob
from the_search.dwarf import Dwarf
from the_search import cuts
import warnings

warnings.filterwarnings("ignore", module='astropy.*')

os.chdir('/Users/runburg/github/gaia_project/dsph_search')

# set number of random coordinates to search through or alternatively amount of candidates to find
num_coords = 1000
num_candids = 100

# import list of previously identified dSph candidates to crossmatch
# open("dsph_candidates.txt")

# append to file if not present

# save new dsph candidate and plots in ...

dwarflist = np.array([
        ['Draco', 260.05972916666667, 57.92121944444444],
        ['LeoI', 152.11716666666666, 12.306500000000002],
        ['LeoB', 168.36720833333334,  22.152805555555553],
        ['UMi', 227.29725, 67.21436111111112],
        ['SclI', 15.039166666666667, -33.70888888888889],
        ['Car', 100.40291666666667, -50.96611111111111],
        ['SxtI', 153.26208333333332, -1.6147222222222224],
        ['FnxI',  39.99708333333333,  -34.44916666666666],
        ['BooI', 210, 14.5],
        ['BooII', 209.5, 12.85],
        ['CetII', 19.475, -17.416666666666668],
        ['ColI', 82.85, -28.033333333333335],
        ['GruII', 331.025,  -46.43333333333333],
        ['CBerI', 186.74583333333334, 23.904166666666665],
        ['Seg1', 151.76333333333335, 16.07361111111111],
        ['Ant2', 143.8868, -36.7673]
    ])
rando_ra = np.array([ 65.47629119236987, 326.48154125527896, 72.7146753092771, 139.1067038659881, 346.25127522798704, 289.1081609393706, 257.0704958874846, 332.56218687772304, 247.96320615330532, 248.8958744426278, 295.0776709675192, 24.216171811340296, 328.3579688732271, 257.77631791021463])
rando_dec = np.array([ 25.257528676684768, 10.474170121853948, -74.84684844097445, -14.851805273629752, 34.65722523571136, -36.555318085426165, 12.406413156130805, -1.2835395220206993, -4.7979387679353325, -20.173861396876788, 48.189809814586845, 34.87627566749339, 27.280421885116603, 63.91444684764842])
rando = np.array(list(zip(rando_ra, rando_dec)))

# dwarflist=dwarflist[8:10]
# rando = []
def create_sample_dwarfs():
    """Sample set of dwarfs for testing."""
    for name, ra, dec in dwarflist:
        Dwarf(ra, dec, name=name, rad_init=0.5)

    for ra, dec in rando:
        Dwarf(ra, dec, rad_init=0.5)


def load_sample_dwarfs(known=True):
        """Load sample set of dwarfs for testing."""
        dwarfs = []
        if known is True:
            for (name, ra, dec) in dwarflist:
                d = Dwarf(ra, dec, name=name)
                for table in glob.glob(f'./candidates/{d.name}/vots/*.vot'):
                    d.load_gaia_table(table)
                yield d
                # dwarfs.append(d)
        else:
            for ra, dec in rando:
                d = Dwarf(ra, dec)
                for table in glob.glob(f'./candidates/{d.name}/vots/*.vot'):
                    d.load_gaia_table(table)
                yield d
                # dwarfs.append(d)

        return dwarfs


def main():
    """Execute cuts."""
    ave = []
    std = []
    cutlist = [5]
    for dwarf in load_sample_dwarfs():
        for cut in cutlist:
            cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True)
            cuts.angular_density_test(dwarf, print_to_stdout=True)


    for dwarf in load_sample_dwarfs(known=False):
        for cut in cutlist:
            cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True)
            cuts.angular_density_test(dwarf, print_to_stdout=True)


if __name__ == "__main__":
    main()
    # create_sample_dwarfs()
    # d = load_sample_dwarfs()
    # dra = Dwarf(260.05972916666667, 57.92121944444444, name='Draco')
    # dra.load_gaia_table('./candidates/Draco/vots/Draco_500.vot')
    # print(dra.gaia_data[-1][-1][[1,2,3]])
