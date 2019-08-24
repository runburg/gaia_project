#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define the Dwarf class.

Author: Jack Runburg
Date: 22-08-2019 14:16


"""
from utils import gaia_search
import os, datetime

class Dwarf:
    """A Dwarf spheroidal galaxy (candidate) object.

    Store relevant dwarf parameters and gaia search information as:
        - name: The given name of a dwarf or defaults to GAIA search id.
        - ra: The right ascension in degrees in ICRS coordinates
        - dec: The declination in degrees in ICRS coordinates
        - gaia_data: The table returned from a GAIA query
        - notes: Note on dwarf or candidacy or crossmatch
    """

    def __init__(self, ra, dec, rad_init=5, name=None):
        """Initialize dwarf object with coordinates."""
        self.ra = ra
        self.dec = dec
        self.gaia_data = []
        self.log = []

        # set name
        if name is not None:
            self.name = name
        else:
            self.name = str(round(ra*100)) + '_' + str(round(dec*100))

        # create home for dwarf candidate data
        os.chdir('./dsph_search')
        os.mkdir(self.name)
        os.mkdir(f'./candidates/{self.name}/plots')
        os.mkdir(f'./candidates/{self.name}/vots')
        self.log.append(f'Created dwarf {self.name} ' + str(datetime.datetime.now()))

        # create GAIA query and store initial data
        self.add_gaia_data(rad_init)


    def add_gaia_data(self, radius):
        """Add gaia search table to Dwarf."""
        job = gaia_search(self.ra, self.dec, self.name, radius)
        self.log.append(f'For radius {radius}; job {job.jobid} stored in {job.outputFile}')
        self.gaia_data.append((radius, job.get_results()))


if __name__ == '__main__':
    rando = Dwarf(12, 20)
    rando.add_gaia_data(0.5)
    # print(rando.gaia_data)
    print('\n'.join(rando.log))
