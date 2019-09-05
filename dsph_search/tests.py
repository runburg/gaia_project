#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit testing for dsph search.

Author: Jack Runburg
Date: 27-08-2019 15:11

"""
import unittest
import os
from astropy.table import Table
from the_search.dwarf import Dwarf

os.chdir('/Users/runburg/github/gaia_project/dsph_search')


class DwarfTests(unittest.TestCase):
    """
    Test the behavior of Dwarf class.

    Ensure no changes to functionality while I fiddle with the class.
    """

    def setUp(self):
        """Create Dwarf object for testing."""
        self.dwarf = Dwarf(0, 10, name='test')

    def test_coords(self):
        """Create object with fed in coordinates."""
        self.assertEqual(self.dwarf.ra, 0)
        self.assertEqual(self.dwarf.dec, 10)
        self.assertEqual(self.dwarf.name, 'test')

    def test_gaia_query(self):
        """Make sure GAIA populated gaia_data dict."""
        self.dwarf.add_gaia_data(0.05)
        self.assertEqual(list(self.dwarf.gaia_data.keys())[-1], 0.05)
        self.assertIsInstance(self.dwarf.gaia_data[0.05], Table)

    def test_rejected(self):
        """Test for cleanup of bad dwarfs."""
        self.dwarf.add_gaia_data(0.005)
        self.dwarf.rejected()
        with self.assertRaises(FileNotFoundError):
            x = open(f'.candidates/{self.dwarf.name}')
        self.assertTrue(os.path.isfile(f'./dead_logs/log_{self.dwarf.name}.log'))


if __name__ == '__main__':
    unittest.main()
