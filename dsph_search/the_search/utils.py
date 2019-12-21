#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for the search.

Author: Jack Runburg
Date: 22-08-2019 14:30


"""

import warnings
from random import random
from time import sleep
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import numpy as np


def gaia_search(ra, dec, name, output_path, radius=0.5, sigma=3, pm_threshold=5, bp_rp_threshold=2, dump_to_file=True):
    """Given coordinates, return gaia cone search around object."""
    warnings.filterwarnings("ignore", module='astropy.*')
    coords = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    job = Gaia.launch_job_async(f"SELECT TOP 500000 \
                                gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec, \
                                gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error, \
                                gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error, \
                                gaia_source.bp_rp, gaia_source.phot_g_mean_mag \
                                FROM gaiadr2.gaia_source \
                                WHERE \
                                CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{coords.ra.degree},{coords.dec.degree},{radius}))=1 AND  (gaiadr2.gaia_source.parallax - gaiadr2.gaia_source.parallax_error * {sigma} <= 0) AND (SQRT(POWER(gaiadr2.gaia_source.pmra, 2) + POWER(gaiadr2.gaia_source.pmdec, 2)) <= {pm_threshold}) AND (gaiadr2.gaia_source.bp_rp <= {bp_rp_threshold})", dump_to_file=dump_to_file, output_file=f'{output_path}/vots/{name}_{round(radius*100)}.vot', verbose=True)

    return job


def get_window_function(spa_dim, pm_dim):
    """Return a window function for a 4d convolution."""
    # ensure int
    spa_dim = int(spa_dim)
    pm_dim = int(pm_dim)

    # create 4d ellipse in histogram space
    a = np.ones((spa_dim*2+1, spa_dim*2+1, pm_dim*2+1, pm_dim*2+1))
    for i in np.arange(-spa_dim, 1+spa_dim):
        for j in np.arange(-spa_dim, spa_dim+1):
            for k in np.arange(-pm_dim, pm_dim+1):
                for l in np.arange(-pm_dim, pm_dim+1):
                    if (i/spa_dim)**2 + (j/spa_dim)**2 + (k/pm_dim)**2 + (l/pm_dim)**2 > 1:
                        a[i+spa_dim, j+spa_dim, k+pm_dim, l+pm_dim] = 0

    return a


def random_cones_outside_galactic_plane(limit=15):
    """Check if in plane and return new coordinates if not."""
    # galactic longitude
    l = random() * 360
    # galactic latitude
    b = (random() - 0.5) * (90 - limit) * 2
    # check that no galactic latitudes are within limit deg of galactic plane

    if b < 0:
        b -= limit
    else:
        b += limit

    c_gal = SkyCoord(l, b, unit='deg', frame='galactic')
    icrs_coords = (c_gal.icrs.ra.value, c_gal.icrs.dec.value)

    return icrs_coords


def gaia_region_search(ra, dec, outfile, radius=10, sigma=3, pm_threshold=5, bp_rp_threshold=2, limit=15,  dump_to_file=True):
    """Given coordinates, search gaia around a region and populate cones within that region."""
    warnings.filterwarnings("ignore", module='astropy.*')
    coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
    job = Gaia.launch_job_async(f"SELECT TOP 10000000 \
                                gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec, \
                                gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error, \
                                gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error, \
                                gaia_source.bp_rp, gaia_source.phot_g_mean_mag \
                                FROM gaiadr2.gaia_source \
                                WHERE \
                                CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',{coords.ra.degree},{coords.dec.degree},{2*radius},{2*radius}))=1 AND  (gaiadr2.gaia_source.parallax - gaiadr2.gaia_source.parallax_error * {sigma} <= 0) AND (SQRT(POWER(gaiadr2.gaia_source.pmra, 2) + POWER(gaiadr2.gaia_source.pmdec, 2)) <= {pm_threshold}) AND (gaiadr2.gaia_source.bp_rp <= {bp_rp_threshold})", dump_to_file=dump_to_file, output_file=outfile, verbose=True)

    return job


def azimuthal_equidistant_coordinates(gaia_table, region_ra, region_dec):
    """Return cartesian coordinates from GAIA table using azimuthal equidistant projection."""
    # http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html

    ra_rad = np.deg2rad(region_ra)
    dec_rad = np.deg2rad(region_dec)

    ra_gaia_rad = np.deg2rad(gaia_table['ra'])
    dec_gaia_rad = np.deg2rad(gaia_table['dec'])

    c = np.arccos(np.sin(dec_rad)*np.sin(dec_gaia_rad) + np.cos(dec_rad)*np.cos(dec_gaia_rad)*np.cos(ra_gaia_rad-ra_rad))

    k_prime = c / np.sin(c)

    x = k_prime * np.cos(dec_gaia_rad) * np.sin(ra_gaia_rad - ra_rad)
    y = k_prime * (np.cos(dec_rad)*np.sin(dec_gaia_rad) - np.sin(dec_rad)*np.cos(dec_gaia_rad)*np.cos(ra_gaia_rad-ra_rad))

    return x, y


def inverse_azimuthal_equidistant_coordinates(x, y, region_ra, region_dec):
    """Given (x, y) positions from AEP, return (ra, dec)."""
    # http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
    ra_rad = np.deg2rad(region_ra)
    dec_rad = np.deg2rad(region_dec)

    c = np.sqrt(x**2 + y**2)

    phi = np.arcsin(np.cos(c)*np.sin(dec_rad) + y/c * np.sin(c)*np.cos(dec_rad))
    lamb = ra_rad + np.arctan2(x*np.sin(c), (c*np.cos(dec_rad)*np.cos(c) - y*np.sin(dec_rad)*np.sin(c)))

    return np.rad2deg(lamb), np.rad2deg(phi)


def get_cone_in_region(ra, dec, region_radius, max_radius=1, limit=15, num_cones=10000):
    """Generate cones given a circular region of region_radius that won't bleed out of the region."""
    sample_region_size = int(np.sqrt(num_cones))
    sample_range_locations = np.linspace(-region_radius+max_radius, region_radius+max_radius, num=sample_region_size)

    x_center = np.sin(np.deg2rad(dec))*np.cos(np.deg2rad(ra))
    y_center = np.sin(np.deg2rad(dec))*np.sin(np.deg2rad(ra))
    z_center = np.cos(np.deg2rad(dec))

    x_locations = x_center + sample_range_locations
    y_locations = y_center + sample_range_locations
    z_locations = np.sqrt((max_radius-region_radius)**2-x_locations**2-y_locations**2)

    theta = np.rad2deg(np.arccos(z_locations/np.sqrt(x_locations**2+y_locations**2+z_locations**2)))
    phi = np.rad2deg(np.arctan2(y_locations, x_locations))

    print([(ph, th) for ph, th in zip(theta, phi)])
    # for point in range(num_cones):
    #     if point % 1000000 == 0:
    #         print(f"At cone {point}")
    #     point += 0.5
    #
    #     # equally spaced coordinates in degrees
    #     theta = 180/np.pi * (np.arccos(1 - 2 * point / num_cones) - np.pi/2)
    #     phi = (180 * (1 + 5**0.5) * point) % 360
    #
    #     # if out of declination range, continue
    #     if abs(dec - theta) > region_radius - max_radius:
    #         continue
    #     else:
    #         ang_dist = angular_distance(ra, dec, phi, theta)
    #         # if outside region, continue
    #         if ang_dist > (region_radius - max_radius) * np.pi/180:
    #             continue
    #         elif outside_of_galactic_plane(phi, theta) is True:
    #             # if outside of galactic plane
    #             yield (phi, theta)


def outside_of_galactic_plane(ra, dec, limit=15):
    """Check that coordinates are outside (up to limit) the galactic plane."""
    c_icrs = SkyCoord(ra, dec, unit='deg', frame='icrs')
    return np.abs(c_icrs.galactic.b.value) > limit


def angular_distance(ra, dec, ra_cone, dec_cone):
    """For two sets of coordinates, find angular_distance between them in radians."""
    ra_diff = ra - ra_cone
    # for i, dif in enumerate(ra_diff):
    if ra_diff < 0:
        ra_diff += 360
    # using vincenty formula from https://en.wikipedia.org/wiki/Great-circle_distance
    ra_diff_rad = abs(np.deg2rad(ra_diff))
    dec_rad = np.deg2rad(dec)
    dec_cone_rad = np.deg2rad(dec_cone)
    # # ang_dist = np.arctan(np.sqrt((np.cos(dec_cone_rad) * np.sin(ra_diff_rad))**2 + (np.cos(dec_rad)*np.sin(dec_cone_rad)-np.sin(dec_rad)*np.cos(dec_cone_rad)*np.cos(ra_diff_rad))**2) /(np.sin(dec_rad)*np.sin(dec_cone_rad)+np.cos(dec_rad)*np.cos(dec_cone_rad)*np.cos(ra_diff_rad)))
    ang_dist = np.arccos(np.sin(dec_rad)*np.sin(dec_cone_rad) + np.cos(dec_rad)*np.cos(dec_cone_rad)*np.cos(ra_diff_rad))

    return ang_dist


def unmask(data):
    """Return an unmasked table for cuts."""
    data = data[[~obj.mask for obj in data]]

    return data


def try_until(func, max_tries=6, sleep_time=30):
    """Try to fetch GAIA data max_tries times."""
    for _ in range(0, max_tries):
        try:
            return func()
        except:
            sleep(sleep_time)
    raise GaiaResultsNotReturnedError()


class GaiaResultsNotReturnedError(Exception):
    """Error for when GAIA is taking too long to return data."""

    pass


def fibonnaci_sphere(num_points, limit=16, point_start=0, point_end=None):
    """Return a coordinate on a Fibonnaci sphere."""
    if point_end is None:
        point_end = num_points
    for point in range(int(point_start), int(point_end)):
        if point % 100 == 0:
            print(point)
        point += 0.5

        # equally spaced coordinates
        theta = 180/np.pi * (np.arccos(1 - 2 * point / num_points) - np.pi/2)
        phi = 180 * (1 + 5**0.5) * point

        if abs(theta) > limit:
            c_gal = SkyCoord(phi, theta, unit='deg', frame='galactic')
            icrs_coords = (c_gal.icrs.ra.value, c_gal.icrs.dec.value)

            # lmc_coords = SkyCoord(80.89, -69.76, unit='deg')
            # lmc = CircleSkyRegion(lmc_coords, Angle(5.5, 'deg'))
            # if icrs_coords not in lmc:
            yield icrs_coords


if __name__ == '__main__':
    # gaia_region_search(90, 90)
    # for _ in range(20):
    #     print('{}, {}'.format(*random_cones_outside_galactic_plane()))
    get_cone_in_region(10, 20, 5, num_cones=20)
