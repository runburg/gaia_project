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


def plot_setup():
    """Set mpl parameters for beautification."""
    return None


def parallax_histograms(input_file, input_pm_file, output_file, titles=None, cut=100):
    """Generate the temperature maps for pms of input objects."""
    # get titles for each subplot and dwarf proper motions
    dwarf_specs = ascii.read(input_file, format='ecsv')
    if titles is None:
        titles = [f'RA: {np.round(ra, 2)}, DEC: {np.round(dec, 2)}' for (ra, dec) in zip(dwarf_specs['RA'], dwarf_specs['DEC'])]
        dwarf_pmra, dwarf_pmdec = None, None
    else:
        titles = dwarf_specs['MAIN_ID']
        dwarf_pmra = dwarf_specs['PMRA']
        dwarf_pmdec = dwarf_specs['PMDEC']

    # load stellar pm values
    with np.load(input_pm_file, 'rb', allow_pickle=True) as infile:
        # vels = infile['vels']
        pmra = infile['pmra']
        pmdec = infile['pmdec']
        parallax = infile['parallax']
        ra = infile['ra']
        dec = infile['dec']

    # setup plot
    plt.close('all')
    mpl.rcParams['axes.labelsize'] = 'x-large'
    mpl.rcParams['ytick.labelsize'] = 'x-small'
    mpl.rcParams['xtick.labelsize'] = 'x-small'

    # set fig size and shape
    rows = 3
    cols = 3
    figsize = (3 * cols, 2 + 2 * rows)
    fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)

    # plot each dwarf in separate subplots
    for i in range(0, rows):
        # set up symmetrical bins about axis
        bins = np.logspace(-2, 0, num=30)

        # add quiver plot
        cut_data = cut_on_parallax(ra[i], dec[i], pmra[i], pmdec[i], cut=cut)
        pm_mag = np.hypot(cut_data[2], cut_data[3])
        arrows = ax[i, 0].quiver(*cut_data, color=cm.gnuplot(pm_mag), clim=(0, pm_mag.max()), units='x', pivot='tail', width=0.005)
        colorbar_for_subplot(fig, ax[i, 0], arrows)
        ax[i, 0].set_facecolor('xkcd:black')
        ax[i, 0].set_title(titles[i].strip('NAME ')+' quiver')

        # generate histogram for each cone
        n, bins, patches = ax[i, 1].hist(parallax[i], bins=bins)
        ax[i, 1].set_title(titles[i].strip('NAME ')+' parallax histogram')
        ax[i, 1].set_xscale('log')

        # add pm 2d histograms
        bins = np.linspace(-15, 15, num=180)
        counts, xedges, yedges, im = ax[i, 2].hist2d(pmra[i], pmdec[i], bins=(bins, bins), vmin=0, cmap='gnuplot')
        if dwarf_pmra is not None:
            ax[i, 2].plot(dwarf_pmra[i], dwarf_pmdec[i], marker='x', color='xkcd:cyan', alpha=0.5)

        ax[i, 2].set_title(titles[i].strip('NAME '))
        colorbar_for_subplot(fig, ax[i, 2], im)

    plt.tight_layout()

    # add a universal colorbar, change cmap in hist2d above
    # fig.colorbar(im, ax=axs.ravel().tolist())

    fig.savefig(output_file, bbox_inches='tight')


def make_pm_maps(input_file, input_pm_file, output_file, num_cones, num_bins=80, titles=None, mincount=0, maxcount=40):
    """Generate the temperature maps for pms of input objects."""
    # get titles for each subplot and dwarf proper motions
    dwarf_specs = ascii.read(input_file, format='ecsv')
    print(dwarf_specs['RA'])
    print(dwarf_specs['DEC'])
    if titles is None:
        titles = [str('RA: ' + str(np.round(ra, 2)) + ', DEC: ' + str(np.round(dec, 2))) for (ra, dec) in zip(dwarf_specs['RA'], dwarf_specs['DEC'])]
        dwarf_pmra, dwarf_pmdec = None, None
    else:
        titles = dwarf_specs['MAIN_ID']
        dwarf_pmra = dwarf_specs['PMRA']
        dwarf_pmdec = dwarf_specs['PMDEC']

    # load stellar pm values
    with np.load(input_pm_file, 'rb', allow_pickle=True) as infile:
        # vels = infile['vels']
        pmra = infile['pmra']
        pmdec = infile['pmdec']

    # setup plot
    plt.close('all')
    mpl.rcParams['axes.labelsize'] = 'x-large'
    mpl.rcParams['ytick.labelsize'] = 'x-small'
    mpl.rcParams['xtick.labelsize'] = 'x-small'

    # set fig size and shape
    d = len(titles)
    rows = 3
    cols = int(np.ceil(d/rows))
    figsize = (4 + 3 * cols, 2 + 2 * rows)
    fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)
    axs = trim_axes(axs, d)
    max_count = [0, 0]

    # plot each dwarf in separate subplots
    for i, ax in enumerate(fig.axes):
        # set up symmetrical bins about axis
        # bins = np.logspace(-2, 2, num=num_bins//2)
        # bins = np.concatenate((np.flip(np.negative(bins)), bins))
        bins = np.linspace(-15, 15, num=180)

        # generate temp map for each cone
        counts, xedges, yedges, im = ax.hist2d(pmra[i], pmdec[i], bins=(bins, bins), vmin=0, cmap='gnuplot')
        if counts.max() > max_count[1]:
            max_count = [i, counts.max()]

        # plot the proper motion of each dwarf from simbad
        if dwarf_pmra is not None:
            ax.plot(dwarf_pmra[i], dwarf_pmdec[i], marker='x', color='xkcd:cyan', alpha=0.5)

        ax.set_title(titles[i].strip('NAME '))
        colorbar_for_subplot(fig, ax, im)
        # ax.set_xscale('symlog', linthreshx=0.1)
        # ax.set_yscale('symlog', linthreshy=0.1)

    # make labels across all subplots
    universal_plot_labels(fig, r"Proper motion, right ascension [mas/yr]", r"Proper motion, declination [mas/yr]")

    plt.tight_layout()

    # add a universal colorbar, change cmap in hist2d above
    # fig.colorbar(im, ax=axs.ravel().tolist())

    fig.savefig(output_file, bbox_inches='tight')
    print('max count is: ' + str(max_count))


def cut_on_parallax(xcoord, ycoord, xpm, ypm, cut):
    """Return objects with proper motion magnitude below cut."""
    indices = [np.hypot(xpp, ypp) < cut for xpp, ypp in zip(xpm,ypm)]
    xpm = xpm[indices]
    ypm = ypm[indices]
    xcoord = xcoord[indices]
    ycoord = ycoord[indices]

    return xcoord, ycoord, xpm, ypm


def universal_plot_labels(fig, xlabel, ylabel):
    """Put big labels across all subplots."""
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()


def colorbar_for_subplot(fig, axs, plot):
    """Place a colorbar by each plot."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(axs)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(plot, cax=cax)


def pm_vector_plot(xcoord, ycoord, xpm, ypm, outfile, names=None, cut=10, cmap=cm.gnuplot):
    """Make vector plots of pm data."""
    d = len(xcoord)
    rows = d//4
    if rows==0:
        rows=1
    cols = int(np.ceil(d/rows))
    figsize = (3.5 * cols, 2.75 * rows)
    fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)
    # fig.tight_layout()
    axs = trim_axes(axs, d)

    for i, (ax, name) in enumerate(zip(axs, names)):
        cut_data = cut_on_parallax(xcoord[i], ycoord[i], xpm[i], ypm[i], cut=cut)
        pm_mag = np.hypot(cut_data[2], cut_data[3])

        arrows = ax.quiver(*cut_data, color=cmap(pm_mag), clim=(0, pm_mag.max()), units='x', pivot='tail', width=0.005)
        # qk = ax3.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')

        colorbar_for_subplot(fig, ax, arrows)
        ax.set_facecolor('xkcd:black')
        ax.set_title(name.strip('NAME '))

    universal_plot_labels(fig, r'Right Ascension [$^\circ$]', r'Declination [$^\circ$]')
    fig.savefig(outfile, bbox_inches='tight')


if __name__ == '__main__':
    dwarf_specs = ascii.read('./gaia_data/dwarf_info.ecsv', format='ecsv')
    xcoords = dwarf_specs['RA']
    ycoords = dwarf_specs['DEC']
    titles = dwarf_specs['MAIN_ID']

    # load pm values
    with np.load('./gaia_data/dwarf_vels.npz', allow_pickle=True) as infile:
        # vels = infile['vels']
        pmra = infile['pmra']
        pmdec = infile['pmdec']
        ra = infile['ra']
        dec = infile['dec']

    pm_vector_plot(ra, dec, pmra, pmdec, './plots/dwarf_quivers.pdf', names=titles)
    plt.show()
