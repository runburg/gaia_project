#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cuts for testing dSph candidacy.

Author: Jack Runburg
Date: 22-08-2019 14:31


"""
import numpy as np
from .utils import unmask
from scipy.signal import find_peaks


def cut_on_parallax(table, cut):
    """Return objects with parallax magnitude below cut."""
    # Check if parallax is consistent with zero, indicating faraway objects
    indices = [(par - par_error * cut) < 0 for par, par_error in zip(table['parallax'].data, table['parallax_error'].data)]
    # indices = np.nonzero(indices)[0]

    return table[indices]


def parallax_test(dwarf, all=True, table_index=-1, cuts=[10,1]):
    """Decide if dwarf passes the parallax test."""
    if all is True:
        tables = dwarf.gaia_data
    else:
        tables = [dwarf.gaia_data[table_index]]

    for table in tables:
        radius = table[0]

        pars = []
        for cut in cuts:
            tablex = cut_on_parallax(table[1], cut)
            parallax = tablex['parallax']
            pars.append(len(parallax))
        print(f'{dwarf.name} @ {radius}: {pars}, {pars[-1]/pars[0]}')


def proper_motion_test(dwarf, all=True, table_index=-1,):
    """Decide if dwarf passes proper motion test."""
    bound = 5
    bins = np.linspace(-bound, bound, num=20*bound)

    if all is True:
        tables = dwarf.gaia_data
    else:
        tables = [dwarf.gaia_data[table_index]]

    for table in tables:
        radius = table[0]
        table = table[1]
        histo, xbins, ybins = np.histogram2d(table['pmra'].data, table['pmdec'].data, bins=(bins, bins))
        # import matplotlib.pyplot as plt
        # plt.hist2d(table['pmra'], table['pmdec'], bins=(bins,bins))
        # plt.show()

    histo_sorted = np.argsort(histo.flatten())
    histo_len = len(histo)
    # max_index = histo.argmax()
    # maximum = histo.max()
    total = np.sum(histo)
    test_area = 15
    test_percentage = 0.3

    for index in range(-10,0):
        max_index = histo_sorted[index]
        coords = (max_index//histo_len, max_index % histo_len)


        log_message = f': Proper motion test: area {test_area*2} x {test_area*2}, percent {test_percentage*100}%.'
        subtotal = np.sum(histo[coords[0]-test_area:coords[0]+test_area+1, coords[1]-test_area:coords[1]+test_area+1])

        if subtotal/total >= test_percentage:
            print(f'{dwarf.name}, proper motion test: PASS')
            dwarf.log.append('PASS' + log_message)
            return True

    print(f'{dwarf.name}, proper motion test: FAIL')
    dwarf.log.append('FAIL' + log_message)
    return False

    # print(f'{dwarf.name}: {find_peaks(histo.flatten(), height=int(0.7*histo.max()), width=3)[0]}')
    # print(f'{dwarf.name}: {histo.max()}')
    # print(f'{dwarf.name}: {maximum}, {refind_max}, {histo.flatten()[max_index-2:max_index+3]}')
