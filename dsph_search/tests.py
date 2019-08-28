#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit testing for dsph search.

Author: Jack Runburg
Date: 27-08-2019 15:11

"""
import unittest
import os
from astropy.table import Table

os.chdir('/Users/runburg/github/gaia_project/dsph_search')

from the_search.dwarf import Dwarf


class DwarfTests(unittest.TestCase):

    def setUp(self):
        self.dwarf = Dwarf(0, 10, name='test')


    def test_coords(self):
        """Create object with fed in coordinates."""
        self.assertEqual(self.dwarf.ra, 0)
        self.assertEqual(self.dwarf.dec, 10)
        self.assertEqual(self.dwarf.name, 'test')


    def test_gaia_query(self):
        self.dwarf.add_gaia_data(0.05)
        self.assertEqual(self.dwarf.gaia_data[0][0], 0.05)
        self.assertIsInstance(self.dwarf.gaia_data[0][1], Table)


    def test_rejected(self):
        self.dwarf.add_gaia_data(0.005)
        self.dwarf.rejected()
        with self.assertRaises(FileNotFoundError):
            open(f'.candidates/{self.dwarf.name}')
        self.assertTrue(os.path.isfile(f'./dead_logs/log_{self.dwarf.name}.log'))


if __name__ == '__main__':
    unittest.main()
