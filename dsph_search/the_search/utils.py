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
                                CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',{coords.ra.degree},{coords.dec.degree},{2*radius},{2*radius}))=1 AND  (gaiadr2.gaia_source.parallax - gaiadr2.gaia_source.parallax_error * {sigma} <= 0) AND (gaiadr2.gaia_source.bp_rp <= {bp_rp_threshold})", dump_to_file=dump_to_file, output_file=f'{output_path}/vots/{name}_{round(radius*100)}.vot', verbose=True)

    return job


def convolve_spatial_histo(gaia_table, region_radius, radii):
    """Convolve the spatial histogram of GAIA data with bin sizes given in radii."""
    from astropy import convolution

    # Bin data at finest resolution
    min_radius = min(radii)
    histo, xedges, yedges = np.histogram2d(gaia_table['x'], gaia_table['y'], bins=region_radius//min_radius)
    # print(histo.shape)
    # put bins in degrees
    xedges *= 180/np.pi
    yedges *= 180/np.pi

    # Set bins for plotting
    X, Y = np.meshgrid(xedges, yedges)
    histo_mask = np.less(X[:-1, :-1]**2 + Y[:-1, :-1]**2, region_radius**2)

    # Convolve the histogram with different size tophats
    convolved_data = []
    for radius in radii:
        convolution_kernel = convolution.Tophat2DKernel(radius//min_radius)
        histo_mask = np.less(X[:-1, :-1]**2 + Y[:-1, :-1]**2, region_radius**2)
        convolved_array = np.multiply(convolution.convolve(histo, convolution_kernel), histo_mask)

        # All convolved data is stored here
        convolved_data.append((radius, convolved_array))

        print(f"finished {radius}")

    return convolved_data, xedges, yedges, X, Y, histo, histo_mask


def convolve_pm_histo(gaia_table, region_radius, radii):
    """Convolve the pm histogram of GAIA data with bin sizes given in radii."""
    from astropy import convolution

    # Bin data at finest resolution
    min_radius = min(radii)
    histo, xedges, yedges = np.histogram2d(gaia_table['pmra'], gaia_table['pmdec'], bins=5//min_radius)

    # Set bins for plotting
    X, Y = np.meshgrid(xedges, yedges)
    histo_mask = np.less(X[:-1, :-1]**2 + Y[:-1, :-1]**2, region_radius**2)

    # Convolve the histogram with different size tophats
    convolved_data = []
    for radius in radii:
        convolution_kernel = convolution.Tophat2DKernel(radius//min_radius)
        histo_mask = np.less(X[:-1, :-1]**2 + Y[:-1, :-1]**2, region_radius**2)
        convolved_array = np.multiply(convolution.convolve(histo, convolution_kernel), histo_mask)

        # All convolved data is stored here
        convolved_data.append((radius, convolved_array))

        print(f"finished {radius}")

    return convolved_data, xedges, yedges, X, Y, histo, histo_mask


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


def inverse_azimuthal_equidistant_coordinates(x, y, ra_rad, dec_rad):
    """Given (x, y) positions from AEP, return (ra, dec) in deg."""
    # http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
    c = np.sqrt(x**2 + y**2)

    phi = np.arcsin(np.cos(c)*np.sin(dec_rad) + y/c * np.sin(c)*np.cos(dec_rad))
    if dec_rad == np.pi/2:
        lamb = ra_rad + np.arctan2(-x, y)
    elif dec_rad == -np.pi/2:
        lamb = ra_rad + np.arctan2(x, y)
    else:
        lamb = ra_rad + np.arctan2(x*np.sin(c), (c*np.cos(dec_rad)*np.cos(c) - y*np.sin(dec_rad)*np.sin(c)))

    return np.rad2deg(lamb), np.rad2deg(phi)


def generate_full_sky_cones(cone_radius, galactic_plane=15, hemi='north', out_to_file=True, output_directory='./region_list/'):
    """Generate full sky coverage of candidate cones."""
    angle = 90 - galactic_plane
    deg_values = np.arange(-angle, angle, cone_radius)
    x, y = np.meshgrid(deg_values, deg_values)

    x = np.concatenate((x.flatten(), (x + cone_radius/2)[:, :-1].flatten()))
    y = np.concatenate((y.flatten(), (y + cone_radius/2)[:-1, :].flatten()))

    inside_of_circle = np.less(x**2+y**2, angle**2)
    x = x[inside_of_circle]
    y = y[inside_of_circle]

    # NORTHERN HEMISPHERE
    ra, dec, ra2, dec2 = [], [], [], []
    if hemi == 'north' or hemi == 'both':
        ra, dec = inverse_azimuthal_equidistant_coordinates(np.deg2rad(x), np.deg2rad(y), 0, np.pi/2)
    # SOUTHERN HEMISPHERE
    if hemi == 'south' or hemi == 'both':
        ra2, dec2 = inverse_azimuthal_equidistant_coordinates(np.deg2rad(x), np.deg2rad(y), 0, -np.pi/2)

    ra_gal = np.concatenate((ra, ra2))
    dec_gal = np.concatenate((dec, dec2))

    ra, dec = galactic_to_icrs(ra_gal, dec_gal)

    if out_to_file is True:
        candidate_per_file = 250
        for i in range(1, len(ra)//candidate_per_file+1):
            with open(output_directory + f"region{i}.txt", 'w') as outfile:
                outfile.write("# candidates for full sky search of dsph\n")
                np.savetxt(outfile, np.array([ra[candidate_per_file*(i-1):candidate_per_file*i], dec[candidate_per_file*(i-1):candidate_per_file*i]]).T, delimiter=" ", comments='#')
    else:
        return ra, dec


def galactic_to_icrs(ra_gal, dec_gal):
    """Return galactic coordinates."""
    from astropy.coordinates import SkyCoord

    coords = SkyCoord(ra_gal, dec_gal, unit='deg', frame='galactic')
    return coords.icrs.ra, coords.icrs.dec


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
    # get_cone_in_region(10, 20, 5, num_cones=20)
    # print(inverse_azimuthal_equidistant_coordinates(np.array([0]), np.array([0.0001]), 0.001, -np.pi/2))
    import matplotlib.pyplot as plt
    from astropy import units as u
    import astropy.coordinates as coord

    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, projection="hammer")
    ra, dec = generate_full_sky_cones(3.16, out_to_file=False, hemi='both')
    ra = coord.Angle(ra)
    ra = ra.wrap_at(180*u.deg)
    dec = coord.Angle(dec)
    ax.scatter(ra.radian, dec.radian, color='xkcd:light grey blue')
    fig.savefig("allsampleconesplot.pdf")
