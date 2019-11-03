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
from the_search.utils import random_cones_outside_galactic_plane, fibonnaci_sphere, get_cone_in_region, gaia_region_search

warnings.filterwarnings("ignore")


def create_sample_dwarfs():
    """Sample set of dwarfs for testing."""
    for ra, dec, name in dwarflist:
        Dwarf(ra, dec, name=name, rad_init=0.5)

    for ra, dec in rando:
        Dwarf(ra, dec, rad_init=0.5)


def load_sample_dwarfs(dwarflist, known=True, path=None):
    """Load sample set of dwarfs for testing."""
    dwarfs = []
    if path is None:
        path = 'candidates'

    if known is True:
        for (ra, dec, name) in dwarflist:
            d = Dwarf(ra, dec, name=name, path=path)
            for table in glob.glob(f'{path}/{d.name}/vots/*.vot'):
                d.load_gaia_table(table)
            yield d
            # dwarfs.append(d)
    else:
        for ra, dec in dwarflist:
            d = Dwarf(ra, dec, path=path)
            for table in glob.glob(f'{path}/{d.name}/vots/*.vot'):
                d.load_gaia_table(table)
            yield d
            # dwarfs.append(d)

    return dwarfs


def look_at_tuned_parameter_values(plot=False):
    """Check params from tuning."""
    dwarflist = []
    rando = []
    for dwa in np.loadtxt('the_search/tuning/tuning_known_dwarfs.txt', dtype=str, delimiter=','):
        dwarflist.append([dwa[1].astype(np.float), dwa[2].astype(np.float), dwa[0]])
    for ran in np.loadtxt('the_search/tuning/tuning_random.txt', delimiter=','):
        rando.append([ran[0].astype(np.float), ran[1].astype(np.float)])

    # Execute cuts.
    cutlist = [100]
    set_of_dwarfs = [dwarflist, rando]
    for list_of_dwarfs, known, label in zip(set_of_dwarfs, [True, False], ['Known', 'Random']):
        dwarfpass = 0
        for dwarf in load_sample_dwarfs(list_of_dwarfs, known=known, path='the_search/tuning/test_candidates'):
            for cut in cutlist:
                pass1 = cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
                pass2 = cuts.angular_density_test(dwarf, print_to_stdout=True, **params)
                dwarfpass += pass1 & pass2

            if plot is True:
                dwarf.accepted(plot)

        print(f'{label} dwarf pass rate \t{dwarfpass}/{len(list_of_dwarfs)}')

    # randompass = 0
    # for dwarf in load_sample_dwarfs(rando, known=, path='the_search/tuning/test_candidates'):
    #     for cut in cutlist:
    #         pass1 = cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
    #         pass2 = cuts.angular_density_test(dwarf, print_to_stdout=True, **params)
    #         randompass += pass1 & pass2
    #
    #     if plot is True:
    #         dwarf.accepted(plot)

    # print(f'Random pass rate \t{randompass}/{len(rando)}')
    print(params)


def write_candidate_coords():
    """Write out positions of dwarf candidates."""
    with open('candidate_coords.txt', 'w') as outfile:
        for file in glob.glob('./candidates/*'):
            ra, _, dec = file.rpartition('/')[-1].partition('_')
            outfile.write(str(round(float(ra)/100, 2)) + ' ' + str(round(float(dec)/100, 2)) + '\n')


def new_main():
    """Search region of sky."""
    region_ra, region_dec = 170, 40
    region_radius = 10
    radii = [(1.5, 0.5), (1.0, 0.1)]

    job = gaia_region_search(region_ra, region_dec)
    gaia_table = job.get_results()

    for coords in get_cone_in_region(region_ra, region_dec, region_radius, num_cones=10000):
        dwa = Dwarf(*coords)

        for large, small in radii:
            dwa.search_loaded_gaia_table(large, gaia_table)
            dwa.search_loaded_gaia_table(small, gaia_table)
        cuts.poisson_overdensity_test(dwa, gaia_table, radii)

        message = ''
        for test, test_name in zip(dwa.tests, ['poisson overdensity test']):
            if test is False:
                message += test_name + 'FAIL'
            else:
                message += test_name + 'PASS'
        if all(dwa.tests):
            dwa.accepted(plot=False, output=False, summary=message, log=True, verbose=False)
        else:
            dwa.rejected(summary=message, log=False)


def main(num_cones, point_start, point_end, plot=False):
    """Run through num_cones to look for candidates."""
    # for _ in range(num_cones):
    #     dwa = Dwarf(*random_cones_outside_galactic_plane())
    for coords in fibonnaci_sphere(num_points=num_cones, point_start=point_start, point_end=point_end):
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
            dwa.accepted(plot=plot, output=False, summary=message, log=False, verbose=False)
        else:
            dwa.rejected(summary=message, log=False)


# params = {'test_area': 10, 'test_percentage': 0.179376451145657, 'num_maxima': 8, 'density_tolerance': 1.362830538392538}

# params = {'test_area': 14, 'test_percentage': 0.32151337896836803, 'num_maxima': 8, 'density_tolerance': 1.269830538392538}

params = {'test_area': 18, 'test_percentage': 0.4547369094279682, 'num_maxima': 8, 'density_tolerance': 1.261830538392538}
# params = {'test_area': 42, 'test_percentage': 0.3380960890954652, 'num_maxima': 8, 'density_tolerance': 1.239830538392538}

if __name__ == "__main__":
    # main(num_cones=10000, point_start=0, point_end=None)
    # write_candidate_coords()
    # create_sample_dwarfs()
    # d = load_sample_dwarfs()
    look_at_tuned_parameter_values()
    # dra = Dwarf(260.05972916666667, 57.92121944444444, name='Draco')
    # dra.load_gaia_table('./candidates/Draco/vots/Draco_500.vot')
    # print(dra.gaia_data[-1][-1][[1,2,3]])
