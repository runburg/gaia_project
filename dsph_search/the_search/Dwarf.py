#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define the Dwarf class.

Author: Jack Runburg
Date: 22-08-2019 14:16


"""
from utils import gaia_search

class Dwarf:
    """A Dwarf spheroidal galaxy (candidate) object.

    Store relevant dwarf parameters and gaia search information as:
        - title: The given title of a dwarf or defaults to GAIA search id.
        - ra: The right ascension in degrees in ICRS coordinates
        - dec: The declination in degrees in ICRS coordinates
        - gaia_data: The table returned from a GAIA query
        - notes: Note on dwarf or candidacy or crossmatch
    """

    def __init__(self, ra, dec, tit=None):
        """Initialize dwarf object with coordinates."""
        self.ra = ra
        self.dec = dec
        self.title = tit
        self._title = self.title
        self.gaia_data = []
        self.notes = ""

    def set_title(self, value):
        if self.title is None:
            self.title = value

    def add_gaia_data(self, radius):
        self.gaia_data = gaia_search(self.ra, self.dec, radius)
        print(self.gaia_data)


if __name__ == '__main__':
    rando = Dwarf(12, 20)
    rando.add_gaia_data(0.5)
    print(rando.gaia_data.get_results())
