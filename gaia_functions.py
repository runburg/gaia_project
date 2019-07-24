#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for GAIA project.

Author: Jack Runburg
Date: 22-06-2019 13:12


"""

import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm


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


def gaia_search(simbad_table):
    """Given a table from Simbad, return gaia cone search around object."""
    from astroquery.gaia import Gaia

    r = []
    for object in simbad_table:
        coords = coord.SkyCoord(ra=object['RA'], dec=object['DEC'], unit=(u.degree, u.degree), frame='icrs')
        # radius = u.Quantity(object['GALDIM_MAJAXIS']/2, u.arcmin).to(u.degree)
        radius = 0.5 * u.degree
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
    axes = axes.flat
    for ax in axes[N:]:
        ax.remove()
    return axes[:N]


def plot_setup():
    """Set mpl parameters for beautification."""
    return None


def make_pm_maps(input_file, input_pm_file, output_file, num_cones, num_bins=80, titles=None, mincount=0, maxcount=40):
    """Generate the temperature maps for pms of input objects."""
    # get titles for each subplot
    dwarf_specs = ascii.read(input_file, format='ecsv')
    if titles is None:
        titles = [str('RA: ' + str(np.round(ra, 2)) + ', DEC: ' + str(np.round(dec, 2))) for (ra, dec) in zip(dwarf_specs['RA'], dwarf_specs['DEC'])]
    else:
        titles = dwarf_specs['MAIN_ID']

    # load pm values
    with np.load(input_pm_file, 'rb', allow_pickle=True) as infile:
        # vels = infile['vels']
        pmra = infile['pmra']
        pmdec = infile['pmdec']

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

    for i, ax in enumerate(fig.axes):
        # set up symmetrical bins about axis
        # bins = np.logspace(-2, 2, num=num_bins//2)
        # bins = np.concatenate((np.flip(np.negative(bins)), bins))
        bins = np.linspace(-15, 15, num=180)

        # generate temp map for each cone
        counts, xedges, yedges, im = ax.hist2d(pmra[i], pmdec[i], bins=(bins, bins), vmin=mincount, vmax=maxcount, cmap='gnuplot')
        if counts.max() > max_count[1]:
            max_count = [i, counts.max()]
        ax.set_title(titles[i].strip('NAME'))
        # ax.set_xscale('symlog', linthreshx=0.1)
        # ax.set_yscale('symlog', linthreshy=0.1)

    # ## big axes
    # hide tick and tick label of the big axes
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)

    # universal axis labels
    plt.xlabel("Proper motion, right ascension [mas/yr]")
    plt.ylabel("Proper motion, declination [mas/yr]")

    plt.tight_layout()

    # add a universal colorbar, change cmap in hist2d above
    fig.colorbar(im, ax=axs.ravel().tolist())

    fig.savefig(output_file, bbox_inches='tight')
    print('max count is: ' + str(max_count))


def pm_vector_plot(xcoord, ycoord, xpm, ypm, outfile, cmap=cm.gnuplot):
    """Make vector plots of pm data."""
    fig, ax = plt.subplots()
    ax.set_title("pivot='tip'; scales with x view")
    # indices = [np.hypot(xpp, ypp) < 0.5 for xpp, ypp in zip(xpm,ypm)]
    # xpm = xpm[indices]
    # ypm = ypm[indices]
    # xcoord = xcoord[indices]
    # ycoord = ycoord[indices]
    arrows = ax.quiver(xcoord, ycoord, xpm, ypm, color=cmap(np.hypot(xpm, ypm)), clim=(0,np.hypot(xpm, ypm).max()), units='x', pivot='tail', width=0.005)
    # qk = ax3.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    ax.set_facecolor('xkcd:black')
    fig.colorbar(arrows)
    fig.savefig(outfile, bbox_inches='tight')


if __name__ == '__main__':
    dwarf_specs = ascii.read('./dwarf_info.ecsv', format='ecsv')
    xcoords = dwarf_specs['RA']
    ycoords = dwarf_specs['DEC']

    # load pm values
    with np.load('/Users/runburg/github/gaia_project/dwarf_vels.npz', allow_pickle=True) as infile:
        # vels = infile['vels']
        pmra = infile['pmra']
        pmdec = infile['pmdec']
        ra = infile['ra']
        dec = infile['dec']

    pm_vector_plot(ra[0], dec[0], pmra[0], pmdec[0], 'draco_quiver.pdf')
    plt.show()
