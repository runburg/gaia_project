#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define the Dwarf class.

Author: Jack Runburg
Date: 22-08-2019 14:16


"""
from astroquery.utils.tap.model.modelutils import read_results_table_from_file
from astropy.io.votable.tree import Table
import os, datetime, glob
try:
    from .utils import gaia_search
except ModuleNotFoundError:
    from utils import gaia_search

os.chdir('/Users/runburg/github/gaia_project/dsph_search')

class Dwarf:
    """A Dwarf spheroidal galaxy (candidate) object.

    Store relevant dwarf parameters and gaia search information as:
        - ra: right ascension in degrees in ICRS coordinates
        - dec: declination in degrees in ICRS coordinates
        - name: given name of a dwarf or defaults to GAIA search id.
        - rad_init: initial radius of cone search for GAIA
        - gaia_data: table returned from a GAIA query
        - log: notes on dwarf or candidacy or crossmatch
        - tests: list of tests performed
    """

    def __init__(self, ra, dec, name=None, rad_init=0):
        """Initialize dwarf object with coordinates."""
        self.ra = ra
        self.dec = dec
        self.gaia_data = []
        self.log = []
        self.tests = []

        # set name
        if name is not None:
            self.name = name
        else:
            self.name = f'{int(round(ra*100))}_{int(round(dec*100))}'

        # create home for dwarf candidate data
        try:
            os.mkdir(f'./candidates/{self.name}')
            os.mkdir(f'./candidates/{self.name}/plots')
            os.mkdir(f'./candidates/{self.name}/vots')
        except FileExistsError:
            pass
        self.log.append(f'Created dwarf {self.name} ' + str(datetime.datetime.now()))

        # create GAIA query and store initial data
        if not rad_init==0:
            self.add_gaia_data(rad_init)


    def add_gaia_data(self, radius):
        """Add gaia search table to Dwarf."""
        job = gaia_search(self.ra, self.dec, self.name, radius)
        self.log.append(f'For radius {radius}; job {job.jobid} stored in {job.outputFile}')
        self.gaia_data.append((radius, job.get_results()))


    def load_gaia_table(self, table):
        """Load a previously saved GAIA .vot table."""
        radius_ = float(table.rpartition('_')[-1].strip('.vot'))/100
        self.gaia_data.append((radius_, read_results_table_from_file(table, output_format='votable', correct_units=True)))
        self.log.append(f'For radius {radius_}; table loaded from {table}')


    def accepted(self):
        """Celebrate a possible dwarf candidate."""
        self.log.append('\n\nACCEPTED')
        self.log.append('Summary: ' + ", ".join(self.tests))

        with open(f'./candidates/{self.name}/log_{self.name}.txt', 'w') as outfile:
            outfile.write("\n".join(self.log))

        print(f'Dwarf {self.name} ACCEPTED')


    def rejected(self, reason=''):
        """Delete rejected dwarf data."""
        self.log.append('\n\nREJECTED')
        self.log.append('Summary: ' + ", ".join(self.tests))
        self.log.append('Rejection reason: ' + reason)

        with open(f'./dead_logs/log_{self.name}.txt', 'w') as outfile:
            outfile.write("\n".join(self.log))

        print(f'Dwarf {self.name} REJECTED')

        for item in os.walk(f'./candidates/{self.name}', topdown=False):
            for file in item[2]:
                os.remove(item[0]+'/'+file)
            for dir in item[1]:
                os.rmdir(item[0]+'/'+dir)

        os.rmdir(f'./candidates/{self.name}')


if __name__ == '__main__':
    rando = Dwarf(12, 20)
    # rando.load_gaia_table('./candidates/1200_2000/vots/1200_2000_50.vot')
    # print(rando.gaia_data[-1])
    rando.add_gaia_data(0.5)
    print(rando.gaia_data[-1])
    print(rando.log)
    # print('\n'.join(rando.log))
