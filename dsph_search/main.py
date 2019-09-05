#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main script for running the dSph search and controlling output.

Author: Jack Runburg
Date: 22-08-2019 14:03

"""

import os
import glob
import warnings
import numpy as np
from the_search.dwarf import Dwarf
from the_search import cuts

warnings.filterwarnings("ignore", module='astropy.*')

os.chdir('/Users/runburg/github/gaia_project/dsph_search')

# set number of random coordinates to search through or alternatively amount of candidates to find
num_coords = 1000
num_candids = 100

# import list of previously identified dSph candidates to crossmatch
# open("dsph_candidates.txt")

# append to file if not present

# save new dsph candidate and plots in ...


def create_sample_dwarfs():
    """Sample set of dwarfs for testing."""
    for ra, dec, name in dwarflist:
        Dwarf(ra, dec, name=name, rad_init=0.5)

    for ra, dec in rando:
        Dwarf(ra, dec, rad_init=0.5)


def load_sample_dwarfs(known=True):
    """Load sample set of dwarfs for testing."""
    dwarfs = []
    if known is True:
        for (ra, dec, name) in dwarflist:
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
    cutlist = [100]
    # params = {'test_area': 5, 'test_percentage': 0.0685171143004674, 'num_maxima': 8, 'density_tolerance': 1.362830538392538}
    params = {'test_area': 10, 'test_percentage': 0.179376451145657, 'num_maxima': 8, 'density_tolerance': 1.362830538392538}
    for dwarf in load_sample_dwarfs():
        for cut in cutlist:
            cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
            cuts.angular_density_test(dwarf, print_to_stdout=True, **params)

    for dwarf in load_sample_dwarfs(known=False):
        for cut in cutlist:
            cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
            cuts.angular_density_test(dwarf, print_to_stdout=True, **params)


if __name__ == "__main__":
    dwarflist = []
    rando = []
    for dwa in np.loadtxt('./the_search/tuning_known_dwarfs.txt', dtype=str, delimiter=','):
        dwarflist.append([dwa[1].astype(np.float), dwa[2].astype(np.float), dwa[0]])
    for ran in np.loadtxt('./the_search/tuning_random.txt', delimiter=','):
        rando.append([ran[0].astype(np.float), ran[1].astype(np.float)])
    main()
    # create_sample_dwarfs()
    # d = load_sample_dwarfs()
    # dra = Dwarf(260.05972916666667, 57.92121944444444, name='Draco')
    # dra.load_gaia_table('./candidates/Draco/vots/Draco_500.vot')
    # print(dra.gaia_data[-1][-1][[1,2,3]])
