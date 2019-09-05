#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run a simulated annealing for parameter space in dSph search.

Author: Jack Runburg
Date: 31-08-2019 16:38


"""

import glob
import warnings
import numpy as np
from simanneal import Annealer
from dsph_search.the_search import cuts
from dsph_search.the_search.dwarf import Dwarf

warnings.filterwarnings("ignore", module='astropy.*')

# os.chdir('/Users/runburg/github/gaia_project/dsph_search')


def truth_to_power(param_dict, kind='pm'):
    """Test values and return a number."""
    # perform tests on all dwarfs and translate to a number
    known_value = []
    random_value = []

    # for each dwarf that passes both tests, increment known value
    if kind == 'pm':
        for dwarf_test in dwarfs:
            known_value.append(cuts.proper_motion_test(dwarf_test, **param_dict))

        # for each random that passes both tests, increment random_value
        for dwarf_test in randoms:
            random_value.append(cuts.proper_motion_test(dwarf_test, **param_dict))

    if kind == 'density':
        for dwarf_test in dwarfs:
            known_value.append(cuts.angular_density_test(dwarf_test, **param_dict))
        for dwarf_test in randoms:
            random_value.append(cuts.angular_density_test(dwarf_test, **param_dict))

    # return the difference; ideally random_value stays at 0 and known_value==len(dwarfs)
    return sum(random_value) - sum(known_value)


class OptimizeCutParameters(Annealer):
    """Find optimum parameters for cuts."""

    # pass extra data (the distance matrix) into the constructor
    def __init__(self, state, dwarf_list, random_list, kind='pm'):
        """Start optimization."""
        self.kind = kind
        self.dwarfs = dwarf_list
        self.randoms = random_list
        super(OptimizeCutParameters, self).__init__(state)  # important!

    def move(self):
        """Choose new values using random walk."""
        if self.kind == 'pm':
            self.state['test_area'] += np.random.randint(-1, high=2)
            self.state['test_percentage'] += np.random.uniform(low=-0.05, high=0.05)
            # self.state['num_maxima'] += int(np.random.uniform(-2, high=3))
        if self.kind == 'density':
            self.state['density_tolerance'] += np.random.uniform(low=-0.05, high=0.05)
        # no efficiency gain, just proof of concept

    def energy(self):
        """Calculate rate of passing."""
        energy = 0
        energy = self.truth_to_power()
        return energy

    def truth_to_power(self):
        """Test values and return a number."""
        # perform tests on all dwarfs and translate to a number
        known_value = 0
        random_value = 0
        param_dict = self.state

        # for each dwarf that passes both tests, increment known value
        if self.kind == 'pm':
            for dwarf_test in dwarfs:
                known_value += cuts.proper_motion_test(dwarf_test, **param_dict)

            # for each random that passes both tests, increment random_value
            for dwarf_test in randoms:
                random_value += cuts.proper_motion_test(dwarf_test, **param_dict)

        if self.kind == 'density':
            for dwarf_test in dwarfs:
                known_value += cuts.angular_density_test(dwarf_test, **param_dict)
            for dwarf_test in randoms:
                random_value += cuts.angular_density_test(dwarf_test, **param_dict)

        # return the difference; ideally random_value stays at 0 and known_value==len(dwarfs)
        result = random_value - known_value

        return result


if __name__ == '__main__':
    # parameter dict
    params = {'test_area': 10, 'test_percentage': 0.179376451145657, 'num_maxima': 8, 'density_tolerance': 1.362830538392538}
    params_values = params.values()
    params_names = params.keys()
    # passing_random_cones = 4

    dwarfs = []
    randoms = []
    for dwa in np.loadtxt('./the_search/tuning_known_dwarfs.txt', dtype=str, delimiter=','):
        d = Dwarf(dwa[1].astype(np.float), dwa[2].astype(np.float), name=dwa[0])
        for table in glob.glob(f'./the_search/test_candidates/{d.name}/vots/*.vot'):
            d.load_gaia_table(table)
        dwarfs.append(d)
    print('loaded dwarfs')
    for ran in np.loadtxt('./the_search/tuning_random.txt', delimiter=','):
        d = Dwarf(ran[0], ran[1])
        for table in glob.glob(f'./the_search/test_candidates/{d.name}/vots/*.vot'):
            d.load_gaia_table(table)
        randoms.append(d)
    print('loaded randoms')
    print(truth_to_power(params))

    ocp = OptimizeCutParameters(params, dwarfs, randoms, kind='pm')
    ocp.Tmax = 1000.0  # Max (starting) temperature
    ocp.Tmin = 5     # Min (ending) temperature
    ocp.steps = 2000   # Number of iterations
    ocp.updates = 10   # Number of updates (by default an update prints to stdout)
    # auto_time = ocp.auto(minutes=0.1, steps=2000)
    # ocp.set_schedule(auto_time)
    # since our state is just a list, slice is the fastest way to copy
    sta, e = ocp.anneal()
    print(f"min energy {e} with state {sta}")
