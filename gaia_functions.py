#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for GAIA project.

Author: Jack Runburg
Date: 22-06-2019 13:12


"""

import astropy.units as u
import astropy.coordinates as coord
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np


def stringdms_to_deg(list):
    """Convert from dms to deg given a list of the three componenets."""
    dd = []
    for item in list:
        val = 0
        for i, num in enumerate(item.split()):
            val += float(num)/(60**i)
        dd.append(val)
    return dd


def stringhms_to_deg(list):
    """Convert from dms to deg given a list of the three componenets."""
    dd = []
    for item in list:
        val = 0
        for i, num in enumerate(item.split()):
            val += float(num)/(60**i)
        dd.append(val*15)
    return dd


def pm_mag(pmdec, pmra):
    """Calculate simple magnitude estimate of proper motion; delete Nan rows as well."""
    ra_index = ~np.isnan(pmra)
    pmdec = pmdec[ra_index]
    pmra = pmra[ra_index]
    dec_index = ~np.isnan(pmdec)
    pmdec = pmdec[dec_index]
    pmra = pmra[dec_index]

    return np.array([np.sqrt(dec**2 + ra**2) for (dec, ra) in zip(pmdec, pmra)])


def reject_outliers(values1, values2, n=5):
    """Reject the outliers of an array that are n sigmas away from the mean."""
    # values1 = values1[~np.isnan(values1)]
    # values2 = values2[~np.isnan(values2)]
    ind1 = abs(values1 - np.mean(values1)) < n * np.std(values1)
    values1 = values1[ind1]
    values2 = values2[ind1]
    ind2 = abs(values2 - np.mean(values2)) < n * np.std(values2)
    values1 = values1[ind2]
    values2 = values2[ind2]
    return values1, values2


def weighted_median(values):
    """Compute the weighted median and return the +/- sigma values."""
    temp_j = np.asarray(values)
    weightsum = temp_j[:, 0].sum()

    # sort, accumulate, and normalize array
    temp_j = temp_j[temp_j[:, 1].argsort()]
    temp_j[:, 0] = np.cumsum(temp_j, axis=0)[:, 0]
    temp_j = np.divide(temp_j, np.array([weightsum, 1]))

    # find the ave and +\- sigma values
    low_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5-.34, side='left')][1]
    ave_jvals = temp_j[np.searchsorted(temp_j[:, 0], 0.5, side='left')][1]
    up_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5+.34, side='right')][1]

    return np.array([low_jvals, ave_jvals, up_jvals])


def weighted_average(listofvalues):
    """Compute weighted average of list (a,b) with 'a' as the weights and 'b' as the values."""
    sum = 0
    weights = 0
    for [w, v] in listofvalues:
        sum += w*v
        weights += w
    return sum/weights


def gaia_search(simbad_table, radius=0.5):
    """Given a table from Simbad, return gaia cone search around object."""
    from astroquery.gaia import Gaia

    radius = radius * u.degree
    r = []
    for object in simbad_table:
        coords = coord.SkyCoord(ra=object['RA'], dec=object['DEC'], unit=(u.degree, u.degree), frame='icrs')
        # radius = u.Quantity(object['GALDIM_MAJAXIS']/2, u.arcmin).to(u.degree)
        j = Gaia.cone_search_async(coords, radius)
        r.append(j.get_results())
    return r


def random_cones_outside_galactic_plane(num, limit=15):
    """Check if in plane and return new coordinates if not."""
    # generate random coords in galactic coordinate system
    rando = np.random.rand(num, 2)
    rando[:, 0] = rando[:, 0] * 360
    rando[:, 1] = (rando[:, 1] - 0.5) * 180

    for i, dec in enumerate(rando[:, 1]):
        # check that no galactic latitudes are within 5 deg of galactic plane
        if np.abs(dec) < limit:
            rando[:, 1][i] = (dec + limit/2) * np.sign(dec)

    icrs_coords = []
    # convert coords to icrs
    for ra, dec in rando:
        c_gal = coord.SkyCoord(l=ra*u.degree, b=dec*u.degree, frame='galactic')
        icrs_coords.append([c_gal.icrs.ra.value, c_gal.icrs.dec.value])

    ra_icrs, dec_icrs = np.hsplit(np.array(icrs_coords), 2)

    return ra_icrs.ravel(), dec_icrs.ravel()


def trim_axes(axes, N):
    """Trim the axes list to proper length."""
    if N>1:
        axes = axes.flat
        for ax in axes[N:]:
            ax.remove()
        return axes[:N]
    else:
        return [axes]


def plot_setup(rows, cols, d=0, buffer=(0, 0)):
    """Set mpl parameters for beautification."""
    # setup plot
    plt.close('all')
    mpl.rcParams['axes.labelsize'] = 'large'
    mpl.rcParams['ytick.labelsize'] = 'x-small'
    mpl.rcParams['xtick.labelsize'] = 'x-small'
    mpl.rcParams['figure.subplot.wspace'] = buffer[0]
    mpl.rcParams['figure.subplot.hspace'] = buffer[1]

    figsize = (6 * cols + buffer[0], 5.5 * rows + buffer[1])
    fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)
    if d is not 0:
        axs = trim_axes(axs, d)

    return fig, axs


def colorbar_for_subplot(fig, axs, cmap, image):
    """Place a colorbar by each plot."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(axs)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(image, cax=cax)
    # cm.ScalarMappable(norm=None, cmap=cmap),

    return cbar


def universal_plot_labels(fig, xlabel, ylabel):
    """Put big labels across all subplots."""
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()


def cut_on_parallax(ra, dec, pmra, pmdec, parallax, cut):
    """Return objects with parallax magnitude below cut."""
    indices = [par < cut for par in parallax]
    pmra = pmra[indices]
    pmdec = pmdec[indices]
    ra = ra[indices]
    dec = dec[indices]
    parallax = parallax[indices]

    return ra, dec, pmra, pmdec, parallax


def load_dwarf_info(file, titles=None):
    """Load the dwarf info from Simbad files."""
    dwarf_specs = ascii.read(file, format='ecsv')


    if titles is None:
        titles = [f'RA: {np.round(ra, 2)}, DEC: {np.round(dec, 2)}' for (ra, dec) in zip(dwarf_specs['RA'], dwarf_specs['DEC'])]
        dwarf_pmra, dwarf_pmdec = [None]*len(titles), [None]*len(titles)
    else:
        titles = dwarf_specs['MAIN_ID']
        dwarf_pmra = dwarf_specs['PMRA']
        dwarf_pmdec = dwarf_specs['PMDEC']

    return np.array(titles), dwarf_pmra, dwarf_pmdec


def load_gaia_search_info(file):
    """Return stellar parameters from gaia search."""
    with np.load(file, 'rb', allow_pickle=True) as infile:
        # vels = infile['vels']
        pmra = infile['pmra']
        pmdec = infile['pmdec']
        parallax = infile['parallax']
        ra = infile['ra']
        dec = infile['dec']

    return ra, dec, pmra, pmdec, parallax


def parallax_histogram(fig, ax, data, title, cut=None):
    """Create a histogram of parallax values."""
    if cut is not None:
        ra, dec, pmra, pmdec, parallax, = cut_on_parallax(*data, cut)
    else:
        ra, dec, pmra, pmdec, parallax, = data

    # bins = np.logspace(-2, 0, num=30)
    bins = np.linspace(parallax.min(), parallax.max(), num = 30)
    n, bins, patches = ax.hist(parallax, bins=bins)
    ax.set_title(title)
    ax.set_xlabel('Parallax [mas]')
    ax.set_ylabel('Bin counts')
    # ax.set_xscale('log')

    return n, bins, patches


def quiver_plot(fig, ax, data, title, cut=None, colorbar=True, radius=None):
    """Create quiver plot of stellar position and proper motion."""
    if cut is not None:
        ra, dec, pmra, pmdec, parallax, = cut_on_parallax(*data, cut)
    else:
        ra, dec, pmra, pmdec, parallax, = data

    pm_mag = np.hypot(pmra, pmdec)
    pmra = pmra / pm_mag
    pmdec = pmdec / pm_mag

    cmap = cm.viridis
    arrows = ax.quiver(ra, dec, pmra, pmdec, color=cmap(pm_mag), clim=(0, pm_mag.max()), units='xy', pivot='tail', width=0.01, headwidth=2, headlength=3, minlength=0.01)
    cbar = colorbar_for_subplot(fig, ax, cmap, image=arrows)
    ax.set_facecolor('xkcd:black')

    if radius is not None:
        ax.set_title(title + f', density: {np.round(len(ra)/(np.pi * radius**2), 2)}')
        with open('densities.txt', 'a') as outfile:
            outfile.write(f'{title}\t {np.round(len(ra)/(np.pi * radius**2), 2)}\n')
    else:
        ax.set_title(title)

    ax.set_xlabel(r"Right ascension [$^\circ$]")
    ax.set_ylabel(r"Declination [$^\circ$]")
    cbar.ax.set_ylabel("Proper motion magnitude [mas/yr]", rotation=270, labelpad=10)


    return arrows


def pm_histogram(fig, ax, data, title, dwarf_pmra=None, dwarf_pmdec=None, cut=None, colorbar=True):
    """Create 2d histogram of proper motions from GAIA search."""
    if cut is not None:
        ra, dec, pmra, pmdec, parallax, = cut_on_parallax(*data, cut)
    else:
        ra, dec, pmra, pmdec, parallax, = data

    # bin data from gaia in 2d histogram
    bound = 5
    bins = np.linspace(-bound, bound, num=20*bound)
    counts, xedges, yedges, im = ax.hist2d(pmra, pmdec, bins=(bins, bins), vmin=0, cmap='gnuplot')

    # plot pm motion of dwarf from simbad
    if dwarf_pmra is not None:
        ax.plot(dwarf_pmra, dwarf_pmdec, marker='x', markersize=15, color='xkcd:cyan', alpha=0.5)

    ax.set_title(title)
    ax.set_xlabel(r"Right ascension proper motion [mas/yr])")
    ax.set_ylabel(r"Declination proper motion [mas/yr]")

    cbar = colorbar_for_subplot(fig, ax, cm.gnuplot, image=im)
    cbar.ax.set_ylabel("Bin counts", rotation=270, labelpad=10)

    return counts, xedges, yedges, im


def draco_comparison_plot(input_file, input_pm_file, output_file, cuts, dwarf=0, titles=None, radius=0):
    """Generate the temperature maps for pms of input objects."""

    # get titles for each subplot and dwarf proper motions
    titles, dwarf_pmra, dwarf_pmdec = load_dwarf_info(input_file, titles)

    # load stellar pm values
    ra, dec, pmra, pmdec, parallax = load_gaia_search_info(input_pm_file)

    # set fig size and shape
    fig, axs = plot_setup(len(cuts), 3)

    # plot each dwarf in separate subplots
    for j, (title, *data) in enumerate(zip(titles, ra, dec, pmra, pmdec, parallax)):
        if j==dwarf:
            for i, cut in enumerate(cuts):
                # add quiver plot
                arrows = quiver_plot(fig, axs[i, 0], data, f"{title.strip('NAME ')} quiver", cut=cut, radius=radius)

                # generate parallax histogram
                n, bins, patches = parallax_histogram(fig, axs[i, 1], data, f"{title.strip('NAME ')} parallax, cut={cut} mas", cut=cut)

                # add pm 2d histograms
                counts, xedges, yedges, im = pm_histogram(fig, axs[i, 2], data, f"{title.strip('NAME ')} pm histo", dwarf_pmra=dwarf_pmra[j], dwarf_pmdec=dwarf_pmdec[j], cut=cut)

    plt.tight_layout()


    fig.savefig(output_file, bbox_inches='tight')


def comparison_plot(input_file, input_pm_file, output_file, titles=None, cut=None):
    """Generate the temperature maps for pms of input objects."""

    # get titles for each subplot and dwarf proper motions
    titles, dwarf_pmra, dwarf_pmdec = load_dwarf_info(input_file, titles)

    # load stellar pm values
    ra, dec, pmra, pmdec, parallax = load_gaia_search_info(input_pm_file)

    # set fig size and shape
    fig, axs = plot_setup(3, 3)

    # choose how many dwarfs to compare
    counter = range(0, 3)

    # plot each dwarf in separate subplots
    for (i, title, *data) in zip(counter, titles, ra, dec, pmra, pmdec, parallax):
        # add quiver plot
        arrows = quiver_plot(fig, axs[i, 0], data, f"{title.strip('NAME ')} quiver", cut=cut)

        # generate parallax histogram
        n, bins, patches = parallax_histogram(fig, axs[i, 1], data, f"{title.strip('NAME ')} parallax", cut=cut)

        # add pm 2d histograms
        counts, xedges, yedges, im = pm_histogram(fig, axs[i, 2], data, f"{title.strip('NAME ')} pm histo", dwarf_pmra=dwarf_pmra[i], dwarf_pmdec=dwarf_pmdec[i], cut=cut)

    plt.tight_layout()


    fig.savefig(output_file, bbox_inches='tight')


def make_pm_maps(input_file, input_pm_file, output_file, num_cones, num_bins=80, titles=None, mincount=0, maxcount=40, cut=None):
    """Generate the temperature maps for pms of input objects."""
    # get titles for each subplot and dwarf proper motions
    titles, dwarf_pmra, dwarf_pmdec, = load_dwarf_info(input_file, titles)

    # load stellar pm values
    ra, dec, pmra, pmdec, parallax, = load_gaia_search_info(input_pm_file)

    # set fig size and shape
    d = len(titles)
    rows = 3
    cols = int(np.ceil(d/rows))
    fig, axs = plot_setup(rows, cols, d)
    max_count = [0, 0]

    # plot each dwarf in separate subplots
    for ax, title, dwarfpmra, dwarfpmdec, *data in zip(axs, titles, dwarf_pmra, dwarf_pmdec, ra, dec, pmra, pmdec, parallax):
        counts, xedges, yedges, im = pm_histogram(fig, ax, data, f"{title.strip('NAME ')}", dwarf_pmra=dwarfpmra, dwarf_pmdec=dwarfpmdec, cut=cut)

    # make labels across all subplots
    universal_plot_labels(fig, r"Proper motion, right ascension [mas/yr]", r"Proper motion, declination [mas/yr]")

    # add a universal colorbar, change cmap in hist2d above
    # fig.colorbar(im, ax=axs.ravel().tolist())

    fig.savefig(output_file, bbox_inches='tight')


def pm_vector_plot(input_file, input_pm_file, outfile, titles=None, cut=None, cmap=cm.gnuplot, radius=None):
    """Make vector plots of pm data."""
    # load dwarf info from simbad
    titles, dwarf_pmra, dwarf_pmdec, = load_dwarf_info(input_file, titles)

    # load pm values
    ra, dec, pmra, pmdec, parallax, = load_gaia_search_info(input_pm_file)

    # set parameters for subplots
    d = len(ra)
    rows = int(np.ceil(d/4))
    cols = int(np.ceil(d/rows))
    fig, axs = plot_setup(rows, cols, d)

    for i, (ax, title, *data) in enumerate(zip(axs, titles, ra, dec, pmra, pmdec, parallax)):
        arrows = quiver_plot(fig, ax, data, title.strip("NAME "), cut=cut, radius=radius)

    universal_plot_labels(fig, r'Right Ascension [$^\circ$]', r'Declination [$^\circ$]')
    fig.savefig(outfile, bbox_inches='tight')


if __name__ == '__main__':
    for radius in [0.5, 0.1, 0.05]:
        pm_vector_plot(f'./gaia_data/randomcone_info_{radius}.ecsv', f'./gaia_data/randomcone_vels_{radius}.npz', f'./plots/random_quivers_{radius}.pdf', titles=None, radius=radius)
