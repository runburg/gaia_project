#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define the Dwarf class.

Author: Jack Runburg
Date: 22-08-2019 14:16


"""

import os
import datetime
import glob
import warnings
from astroquery.utils.tap.model.modelutils import read_results_table_from_file
from astropy.io.votable.tree import Table
try:
    from the_search.utils import gaia_search
    from the_search.plots import pm_histogram, parallax_histogram, quiver_plot
except ModuleNotFoundError:
    from utils import gaia_search
    from plots import pm_histogram, parallax_histogram, quiver_plot

warnings.filterwarnings("ignore", module='astropy.*')
os.chdir('/Users/runburg/github/gaia_project/dsph_search')


class Dwarf:
    """A Dwarf spheroidal galaxy (candidate) object.

    Store relevant dwarf parameters and gaia search information as:
        - ra: right ascension in degrees in ICRS coordinates
        - dec: declination in degrees in ICRS coordinates
        - name: given name of a dwarf or defaults to GAIA search id.
        - rad_init: initial radius of cone search for GAIA
        - gaia_data: tuple of search radius and table returned from a GAIA query
        - log: notes on dwarf or candidacy or crossmatch
        - tests: list of tests performed
    """

    def __init__(self, ra, dec, name=None, rad_init=0):
        """Initialize dwarf object with coordinates."""
        self.ra = ra
        self.dec = dec
        self.gaia_data = {}
        self.log = []
        self.tests = []

        # set name
        if name is not None:
            self.name = name
        else:
            self.name = f'{int(round(ra*100))}_{int(round(dec*100))}'

        self.path = f'./candidates/{self.name}'
        # create home for dwarf candidate data
        try:
            os.mkdir(self.path)
            os.mkdir(f'{self.path}/plots')
            os.mkdir(f'{self.path}/vots')
        except FileExistsError:
            pass
        self.log.append(f'Created dwarf {self.name} in {self.path}' + str(datetime.datetime.now()))

        # create GAIA query and store initial data
        if rad_init != 0:
            self.add_gaia_data(rad_init)

    def add_gaia_data(self, radius):
        """Add gaia search table to Dwarf."""
        # automatically cut on parallax with sigma=5 and on pm with pm_threshold=5
        job = gaia_search(self.ra, self.dec, self.name, radius)
        self.log.append(f'For radius {radius}; job {job.jobid} stored in {job.outputFile}')
        table = job.get_results()

        self.gaia_data[radius] = table

    def load_gaia_table(self, table):
        """Load a previously saved GAIA .vot table."""
        # if this line breaks, it's fucked
        # apparently very shallow wrap of result = astropy.table.Table.read(file_name, format=output_format)
        data = read_results_table_from_file(table, output_format='votable', correct_units=True)
        radius = float(table.rpartition('_')[-1].strip('.vot')) / 100

        self.gaia_data[radius] = data
        self.log.append(f'For radius {radius}; table loaded from {table}')

    def accepted(self, output=False, summary=''):
        """Celebrate a possible dwarf candidate."""
        self.log.append('\n\nACCEPTED')
        self.log.append('Summary: ' + summary)

        with open(f'{self.path}/log_{self.name}.txt', 'w') as outfile:
            outfile.write("\n".join(self.log))

        parallax_histogram(self)
        pm_histogram(self)
        quiver_plot(self)

        if output is True:
            print(f'Dwarf {self.name} ACCEPTED')

    def rejected(self, output=False, summary=''):
        """Delete rejected dwarf data."""
        self.log.append('\n\nREJECTED')
        self.log.append('Summary: ' + summary)

        with open(f'./dead_logs/log_{self.name}.log', 'w') as outfile:
            outfile.write("\n".join(self.log))

        if output is True:
            print(f'Dwarf {self.name} REJECTED')

        for item in os.walk(self.path, topdown=False):
            for file in item[2]:
                os.remove(item[0] + '/' + file)
            for folder in item[1]:
                os.rmdir(item[0] + '/' + folder)

        os.rmdir(self.path)


if __name__ == '__main__':
    rando = Dwarf(12, 20)
    # rando.load_gaia_table('./candidates/1200_2000/vots/1200_2000_50.vot')
    # print(rando.gaia_data[-1])
    rando.add_gaia_data(0.5)
    print(rando.gaia_data[0.5])
    print(rando.log)
    # print('\n'.join(rando.log))
