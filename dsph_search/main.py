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
from the_search.utils import random_cones_outside_galactic_plane, fibonnaci_sphere

warnings.filterwarnings("ignore", module='astropy.*')


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


def look_at_tuned_parameter_values():
    """Check params from tuning."""
    dwarflist = []
    rando = []
    for dwa in np.loadtxt('./the_search/tuning/tuning_known_dwarfs.txt', dtype=str, delimiter=','):
        dwarflist.append([dwa[1].astype(np.float), dwa[2].astype(np.float), dwa[0]])
    for ran in np.loadtxt('./the_search/tuning/tuning_random.txt', delimiter=','):
        rando.append([ran[0].astype(np.float), ran[1].astype(np.float)])

    # Execute cuts.
    cutlist = [100]
    # params = {'test_area': 5, 'test_percentage': 0.0685171143004674, 'num_maxima': 8, 'density_tolerance': 1.362830538392538}
    for dwarf in load_sample_dwarfs():
        for cut in cutlist:
            cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
            cuts.angular_density_test(dwarf, print_to_stdout=True, **params)

    for dwarf in load_sample_dwarfs(known=False):
        for cut in cutlist:
            cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
            cuts.angular_density_test(dwarf, print_to_stdout=True, **params)


def main(num_cones, point_start, point_end):
    """Run through num_ocones to look for candidates."""
    # for _ in range(num_cones):
    #     dwa = Dwarf(*random_cones_outside_galactic_plane())
    for coords in fibonnaci_sphere(num_cones, point_start, point_end):
        dwa = Dwarf(*coords)
        cuts.proper_motion_test(dwa, **params)
        cuts.angular_density_test(dwa, **params)

        message = ''
        for test, test_name in zip(dwa.tests, [' pm test ', ' ang.den. ']):
            if test is False:
                message += test_name + 'FAIL'
            else:
                message += test_name + 'PASS'
        if all(dwa.tests):
            dwa.accepted(summary=message)
        else:
            dwa.rejected(summary=message)

    with open('candidate_coords.txt', 'a') as outfile:
        for file in glob.glob('./candidates/*'):
            ra, _, dec = file.rpartition('/')[-1].partition('_')
            outfile.write(str(round(float(ra)/100, 2)) + '\t' + str(round(float(dec)/100, 2)) + '\n')


params = {'test_area': 10, 'test_percentage': 0.179376451145657, 'num_maxima': 8, 'density_tolerance': 1.362830538392538}
if __name__ == "__main__":
    main(num_cones=10000, point_start=0, point_end=None)

    # create_sample_dwarfs()
    # d = load_sample_dwarfs()
    # dra = Dwarf(260.05972916666667, 57.92121944444444, name='Draco')
    # dra.load_gaia_table('./candidates/Draco/vots/Draco_500.vot')
    # print(dra.gaia_data[-1][-1][[1,2,3]])
