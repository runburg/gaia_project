#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main script for running the dSph search and controlling output.

Author: Jack Runburg
Date: 22-08-2019 14:03

"""

import glob
import warnings
import numpy as np
import sys
from astropy.table import Table
from the_search.dwarf import Dwarf
from the_search import cuts
from the_search.utils import fibonnaci_sphere, get_cone_in_region, gaia_region_search, outside_of_galactic_plane, azimuthal_equidistant_coordinates, inverse_azimuthal_equidistant_coordinates, get_window_function
from the_search.plots import get_points_of_circle, convolved_histograms, convolved_histograms_1d, new_all_sky

warnings.filterwarnings("ignore")


def create_sample_dwarfs():
    """Sample set of dwarfs for testing."""
    for ra, dec, name in dwarflist:
        Dwarf(ra, dec, name=name, rad_init=0.5)

    for ra, dec in rando:
        Dwarf(ra, dec, rad_init=0.5)


def load_sample_dwarfs(dwarflist, known=True, path=None):
    """Load sample set of dwarfs for testing."""
    dwarfs = []
    if path is None:
        path = 'candidates'

    if known is True:
        for (ra, dec, name) in dwarflist:
            d = Dwarf(ra, dec, name=name, path=path)
            for table in glob.glob(f'{path}/{d.name}/vots/*.vot'):
                d.load_gaia_table(table)
            yield d
            # dwarfs.append(d)
    else:
        for ra, dec in dwarflist:
            d = Dwarf(ra, dec, path=path)
            for table in glob.glob(f'{path}/{d.name}/vots/*.vot'):
                d.load_gaia_table(table)
            yield d
            # dwarfs.append(d)

    return dwarfs


def look_at_tuned_parameter_values(plot=False):
    """Check params from tuning."""
    dwarflist = []
    rando = []
    for dwa in np.loadtxt('the_search/tuning/tuning_known_dwarfs_old_names.txt', dtype=str, delimiter=','):
        dwarflist.append([dwa[1].astype(np.float), dwa[2].astype(np.float), dwa[0]])
    for ran in np.loadtxt('./dsph_search/the_search/tuning/tuning_random.txt', delimiter=','):
        rando.append([ran[0].astype(np.float), ran[1].astype(np.float)])

    # Execute cuts.
    cutlist = [100]
    set_of_dwarfs = [dwarflist, rando]
    for list_of_dwarfs, known, label in zip(set_of_dwarfs, [True, False], ['Known', 'Random']):
        dwarfpass = 0
        for dwarf in load_sample_dwarfs(list_of_dwarfs, known=known, path='./dsph_search/the_search/tuning/test_candidates'):
            for cut in cutlist:
                pass1 = cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
                pass2 = cuts.angular_density_test(dwarf, print_to_stdout=True, **params)
                dwarfpass += pass1 & pass2

            if plot is True:
                dwarf.accepted(plot)

        print(f'{label} dwarf pass rate \t{dwarfpass}/{len(list_of_dwarfs)}')

    # randompass = 0
    # for dwarf in load_sample_dwarfs(rando, known=, path='the_search/tuning/test_candidates'):
    #     for cut in cutlist:
    #         pass1 = cuts.proper_motion_test(dwarf, cut=cut, print_to_stdout=True, **params)
    #         pass2 = cuts.angular_density_test(dwarf, print_to_stdout=True, **params)
    #         randompass += pass1 & pass2
    #
    #     if plot is True:
    #         dwarf.accepted(plot)

    # print(f'Random pass rate \t{randompass}/{len(rando)}')
    print(params)


def write_candidate_coords():
    """Write out positions of dwarf candidates."""
    with open('candidate_coords.txt', 'w') as outfile:
        for file in glob.glob('./candidates/*'):
            ra, _, dec = file.rpartition('/')[-1].partition('_')
            outfile.write(str(round(float(ra)/100, 2)) + ' ' + str(round(float(dec)/100, 2)) + '\n')


def new_main(param_args):
    """Search region of sky."""
    # set parameters of search
    if len(param_args) > 1:
        name = param_args[1]
        region_ra, region_dec, region_radius, num_cones, *radii = [float(arf) for arf in param_args[2:]]
        num_cones = int(num_cones)
    else:
        region_ra, region_dec = 250, 60
        region_radius = 5
        # 15 degree close to galactic plane takes ~60 min
        num_cones = 10000000
        radii = [1.5, 1.0, 0.5]

    minimum_count = 2
    sigma_threshhold = 3

    # standard paths
    infile = f'regions/region_ra{round(region_ra*100)}_dec{round(region_dec*100)}_rad{round(region_radius*100)}.vot'
    outfile = f'region_candidates/region_ra{round(region_ra*100)}_dec{round(region_dec*100)}_rad{round(region_radius*100)}_candidates.txt'

    with open(outfile, 'a') as fil:
        fil.write(f'# successful candidates for region at ({region_ra}, {region_dec}) and radius {region_radius}')

    # first try to find file
    try:
        gaia_table = Table.read(infile, format='votable')
        print(f"Table loaded from: {infile}")
        print(f"Number of objects: {len(gaia_table)}")
    except FileNotFoundError:
        job = gaia_region_search(region_ra, region_dec, outfile=infile, radius=region_radius)
        gaia_table = job.get_results()
        print("Finished querying Gaia")
        gaia_table = gaia_table[outside_of_galactic_plane(gaia_table['ra'], gaia_table['dec'])]
        print("Finished filtering Gaia table")
        gaia_table['x'], gaia_table['y'] = azimuthal_equidistant_coordinates(gaia_table, region_ra, region_dec)
        print("Finished calculating x-y values done")
        print(f"Table dumped to: {infile}")
        print(f"Number of objects: {len(gaia_table)}")
        gaia_table.write(infile, overwrite='True', format='votable')

    # ## NEW STUFF
    from astropy import convolution

    # bin data at finest resolution
    min_radius = min(radii)
    histo, xedges, yedges = np.histogram2d(gaia_table['x'], gaia_table['y'], bins=region_radius//min_radius)
    # print(histo.shape)
    # put bins in degrees
    xedges *= 180/np.pi
    yedges *= 180/np.pi

    # set bins for plotting
    X, Y = np.meshgrid(xedges, yedges)
    histo_mask = np.less(X[:-1, :-1]**2 + Y[:-1, :-1]**2, region_radius**2)

    # convolve the histogram with different size tophats
    convolved_data = []
    for radius in radii:
        convolution_kernel = convolution.Tophat2DKernel(radius//min_radius)
        # not_mask = np.logical_not(histo_mask)
        histo_mask = np.less(X[:-1, :-1]**2 + Y[:-1, :-1]**2, region_radius**2)
        convolved_array = np.multiply(convolution.convolve(histo, convolution_kernel), histo_mask)
        # print(convolved_array.shape)
        convolved_data.append((radius, convolved_array))
        print(f"finished {radius}")

    passing_xy = cuts.histogram_overdensity_test(convolved_data, (xedges, yedges, histo), region_ra, region_dec, outfile, histo_mask, num_sigma=2, repetition=2)

    # plot the convolved data
    convolved_histograms(convolved_data, (X, Y, histo), passingxy=passing_xy, name=name, region_radius=region_radius)
    convolved_histograms_1d(convolved_data, (X, Y, histo), name=name, mask=histo_mask, region_radius=region_radius)

    # ## END NEW STUFF

####
    # for radius in radii:
    #     poisson_sd = np.sqrt(len(gaia_table) * radius**2/region_radius**2)
    #     print(poisson_sd)
    #     histo, xedges, yedges = np.histogram2d(gaia_table['x'], gaia_table['y'], bins=region_radius//radius)
    #     bin_width = (xedges[1]-xedges[0])/2
    #
    #     passing_indices_y, passing_indices_x = np.argwhere(np.logical_and(np.less(poisson_sd*sigma_threshhold, histo), histo > minimum_count)).T
    #     passing_ra, passing_dec = inverse_azimuthal_equidistant_coordinates(xedges[passing_indices_x]+bin_width, yedges[passing_indices_y]+bin_width, region_ra, region_dec)
    #
    #     with open(outfile, 'w') as outfl:
    #         for ra, dec in zip(passing_ra, passing_dec):
    #             outfl.write(f"{ra} {dec}")
    #
    # for radius in radii:
    #     fig, ax = plt.subplots()
    #     ax.hist2d(gaia_table['x'], gaia_table['y'], bins=region_radius//radius)
    #     fig.savefig(f'sculptor_plots/sculptor_histo_{radius}.pdf')
    #
    # fig, ax = plt.subplots()
    # ax.scatter(passing_ra, passing_dec, s=1)
    # ax.set_xlim(left=12.4, right=17.5)
    # ax.set_ylim(bottom=-30.5, top=-36.9)
    # # ax.scatter(*get_points_of_circle(region_ra, region_dec, region_radius).T, s=1)
    # fig.savefig('sculptor_plots/sculptor_passing_coords_001.pdf')

    # for coords in get_cone_in_region(region_ra, region_dec, region_radius, num_cones=num_cones):
    #     # print(f"found coords {coords}")
    #     dwa = Dwarf(*coords)
    #
    #     dwa.search_loaded_gaia_table(radii, gaia_table)
    #
    #     # print("filled tables")
    #     cuts.poisson_overdensity_test(dwa, gaia_table, region_radius)
    #     # print("finished cut")
    #
    #     message = ''
    #     for test, test_name in zip(dwa.tests, ['poisson overdensity test']):
    #         if test is False:
    #             message += test_name + 'FAIL'
    #         else:
    #             message += test_name + 'PASS'
    #     if all(dwa.tests):
    #         dwa.accepted(plot=False, output=False, summary=message, log=False, verbose=False, coord_file_path=outfile)
    #         # print("passed!")
    #     else:
    #         dwa.rejected(summary=message, log=False)
    #        # print("failed")
####

def random_poisson(param_args):
    """Generate random poisson GAIA data for testing."""
    import matplotlib.pyplot as plt
    # # test_cases = [(radius1, radius2) for radius1 in radii for radius2 in radii+[region_radius] if radius1 < radius2]
    # # print(test_cases)

    region_ra, region_dec, region_radius, num_cones, *radii = [float(arf) for arf in param_args[1:]]

    radii = [0.01, 0.005, 0.001, 0.0005, 0.0001]
    region_radius = 0.05
    num_pts = 20000

    sigma_threshhold = 5
    minimum_count = 3
    x = np.random.uniform(-region_radius, region_radius, num_pts)
    y = np.random.uniform(-region_radius, region_radius, num_pts)

    for radius in radii:
        poisson_sd = np.sqrt(num_pts * radius**2/region_radius**2)
        print(poisson_sd)
        histo, xedges, yedges = np.histogram2d(x, y, bins=region_radius//radius)
        bin_width = (xedges[1]-xedges[0])/2
        passing_indices_y, passing_indices_x = np.argwhere((histo > poisson_sd*sigma_threshhold) & (histo > minimum_count)).T
        passing_ra, passing_dec = xedges[passing_indices_x]+bin_width, yedges[passing_indices_y]+bin_width

        # with open(outfile, 'w') as outfl:
        #     for ra, dec in zip(passing_ra, passing_dec):
        #         outfl.write(f"{ra} {dec}")

    for radius in radii:
        fig, ax = plt.subplots()
        ax.hist2d(x, y, bins=region_radius//radius)
        fig.savefig(f'poisson_simulation/poisson_plot_{radius}.pdf')

    fig, ax = plt.subplots()
    ax.scatter(passing_ra, passing_dec, s=1)
    # ax.scatter(*get_points_of_circle(region_ra, region_dec, region_radius).T, s=1)
    fig.savefig('poisson_simulation/poisson_passing_coords_001.pdf')


def main(num_cones=1000, point_start=0, point_end=None, plot=False):
    """Run through num_cones to look for candidates."""
    # for _ in range(num_cones):
    #     dwa = Dwarf(*random_cones_outside_galactic_plane())
    for coords in fibonnaci_sphere(num_points=num_cones, point_start=point_start, point_end=point_end):
        dwa = Dwarf(*coords)
        cuts.proper_motion_test(dwa, **params)
        cuts.angular_density_test(dwa, **params)

        message = ''
        for test, test_name in zip(dwa.tests, [' pm test ', ' ang.den. ']):
            if test is False:
                message += test_name + 'FAIL'
            else:
                message += test_name + 'PASS'
        if all(dwa.tests):
            dwa.accepted(plot=plot, output=False, summary=message, log=False, verbose=False)
        else:
            dwa.rejected(summary=message, log=False)


def main_4d(param_args):
    """Search region of sky."""
    # set parameters of search
    if len(param_args) > 1:
        name = param_args[1]
        region_ra, region_dec, region_radius, num_cones, radii, pm_spacing = param_args[2:]
        num_cones = int(num_cones)
    else:
        region_ra, region_dec = 250, 60
        region_radius = 5
        # 15 degree close to galactic plane takes ~60 min
        num_cones = 10000000
        radii = [1.5, 1.0, 0.5]

    minimum_count = 2
    sigma_threshhold = 3
    pm_threshhold = 5

    # standard paths
    infile = f'regions/region_ra{round(region_ra*100)}_dec{round(region_dec*100)}_rad{round(region_radius*100)}.vot'
    outfile = f'region_candidates/region_ra{round(region_ra*100)}_dec{round(region_dec*100)}_rad{round(region_radius*100)}_candidates.txt'

    with open(outfile, 'a') as fil:
        fil.write(f'# successful candidates for region at ({region_ra}, {region_dec}) and radius {region_radius}')

    # first try to find file
    try:
        gaia_table = Table.read(infile, format='votable')
        print(f"Table loaded from: {infile}")
        print(f"Number of objects: {len(gaia_table)}")
    except FileNotFoundError:
        job = gaia_region_search(region_ra, region_dec, outfile=infile, radius=region_radius)
        gaia_table = job.get_results()
        print("Finished querying Gaia")
        gaia_table = gaia_table[outside_of_galactic_plane(gaia_table['ra'], gaia_table['dec'])]
        print("Finished filtering Gaia table")
        gaia_table['x'], gaia_table['y'] = azimuthal_equidistant_coordinates(gaia_table, region_ra, region_dec)
        print("Finished calculating x-y values done")
        print(f"Table dumped to: {infile}")
        print(f"Number of objects: {len(gaia_table)}")
        gaia_table.write(infile, overwrite='True', format='votable')

    # ## NEW STUFF
    # from astropy import convolution
    from scipy.signal import convolve

    # bin data at finest resolution
    min_radius = min(radii)
    min_pm_spacing = min(pm_spacing)
    print(min_radius, min_pm_spacing)

    bins = [region_radius//min_radius]*2 + [pm_threshhold//min_pm_spacing]*2
    data_4d = np.array([gaia_table['x'], gaia_table['y'], gaia_table['ra'], gaia_table['dec']]).T
    histo, edges = np.histogramdd(data_4d, bins=bins)
    print(histo.shape)

    # put bins in degrees
    edges[0] *= 180/np.pi
    edges[1] *= 180/np.pi

    # set bins for plotting
    X, Y = np.meshgrid(edges[0], edges[1])
    print(X.shape)
    # histo_mask = np.less(X[:-1, :-1]**2 + Y[:-1, :-1]**2, region_radius**2)

    # convolve the histogram with different size tophats
    convolved_data = []
    for radius, pm_space in zip(radii, pm_spacing):
        window_function = get_window_function(radius//min_radius, pm_space//min_pm_spacing)
        convolved_array = convolve(histo, window_function)
        convolved_data.append((radius, convolved_array))
        print(f"finished {radius}")

    # passing_xy = cuts.histogram_overdensity_test(convolved_data, (xedges, yedges, histo), region_ra, region_dec, outfile, histo_mask, num_sigma=2, repetition=2)
    #
    # # plot the convolved data
    # convolved_histograms(convolved_data, (X, Y, histo), passingxy=passing_xy, name=name, region_radius=region_radius)
    # convolved_histograms_1d(convolved_data, (X, Y, histo), name=name, mask=histo_mask, region_radius=region_radius)


# params = {'test_area': 10, 'test_percentage': 0.179376451145657, 'num_maxima': 8, 'density_tolerance': 1.362830538392538}

# params = {'test_area': 14, 'test_percentage': 0.32151337896836803, 'num_maxima': 8, 'density_tolerance': 1.269830538392538}

params = {'test_area': 18, 'test_percentage': 0.4547369094279682, 'num_maxima': 8, 'density_tolerance': 1.261830538392538}
# params = {'test_area': 42, 'test_percentage': 0.3380960890954652, 'num_maxima': 8, 'density_tolerance': 1.239830538392538}

if __name__ == "__main__":
    import time
    start_time = time.time()
    # random_poisson(sys.argv)
    if len(sys.argv) < 2:
        names = ['Crater2', 'Sculptor', 'Draco', 'HydrusI']
        coords = [(177.3100, -18.413), (15.0392, -33.7089), (260.059728, 57.921219), (37.3890, -79.3089)]
        region_radius = 1
        radii = [0.316, 0.1, 0.0316, 0.01]
        pm_spacing = np.array(radii) * 10
        num_cones = 10

        dwarfs = np.loadtxt('./the_search/tuning/tuning_known_dwarfs.txt', delimiter=",", dtype=str)
        # for name, (ra, dec) in zip(names, coords):
        for name, ra, dec in dwarfs[:]:
            print(name)
            ra = float(ra)
            dec = float(dec)
            param_args = [0, name, ra, dec, region_radius, num_cones] + [radii] + [pm_spacing]
            main_4d(param_args)
            print(f"finished with dwarf {name}\n\n\n")

        new_all_sky(region_radius)
    else:
        new_main(sys.argv)

    print("--- %s seconds ---" % (time.time() - start_time))

    # main(num_cones=10000, point_start=0, point_end=None)
    # write_candidate_coords()
    # create_sample_dwarfs()
    # d = load_sample_dwarfs()
    # look_at_tuned_parameter_values()
    # dra = Dwarf(260.05972916666667, 57.92121944444444, name='Draco')
    # dra.load_gaia_table('./candidates/Draco/vots/Draco_500.vot')
    # print(dra.gaia_data[-1][-1][[1,2,3]])
