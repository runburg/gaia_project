#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for the search.

Author: Jack Runburg
Date: 22-08-2019 14:30


"""

import warnings
from random import random
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def gaia_search(ra, dec, name, radius=0.5, sigma=5, pm_threshold=5, dump_to_file=True):
    """Given a table from Simbad, return gaia cone search around object."""
    warnings.filterwarnings("ignore", module='astropy.*')
    # radius = radius * u.degree
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    # radius = u.Quantity(object['GALDIM_MAJAXIS']/2, u.arcmin).to(u.degree)
    job = Gaia.launch_job_async(f"SELECT TOP 50000 \
                                gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec, \
                                gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error, \
                                gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error, \
                                gaia_source.bp_rp \
                                FROM gaiadr2.gaia_source \
                                WHERE \
                                CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{coords.ra.degree},{coords.dec.degree},{radius}))=1 AND  (gaiadr2.gaia_source.parallax - gaiadr2.gaia_source.parallax_error * {sigma} <= 0) AND (SQRT(POWER(gaiadr2.gaia_source.pmra, 2) + POWER(gaiadr2.gaia_source.pmdec, 2)) <= {pm_threshold})", dump_to_file=dump_to_file, output_file=f'./candidates/{name}/vots/{name}_{round(radius*100)}.vot')

    return job


def random_cones_outside_galactic_plane(limit=15):
    """Check if in plane and return new coordinates if not."""
    # galactic longitude
    l = random() * 360
    # galactic latitude
    b = (random() - 0.5) * (90 - limit)
    # check that no galactic latitudes are within limit deg of galactic plane

    if b < 0:
        b -= limit
    else:
        b += limit

    c_gal = SkyCoord(l=l * u.degree, b=b * u.degree, frame='galactic')
    icrs_coords = (c_gal.icrs.ra.value, c_gal.icrs.dec.value)

    return icrs_coords


def unmask(data):
    """Return an unmasked table for cuts."""
    data = data[[~obj.mask for obj in data]]

    return data


def colorbar_for_subplot(fig, axs, cmap, image):
    """Place a colorbar by each plot."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(axs)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(image, cax=cax)
    # cm.ScalarMappable(norm=None, cmap=cmap),

    return cbar


def parallax_histogram(dwarf):
    """Create a histogram of parallax values."""
    fig, axs = plt.subplots(ncols=2, nrows=len(dwarf.gaia_data)//2 + len(dwarf.gaia_data) % 2)

    for ax, table in zip(axs, dwarf.gaia_data):
        parallax = dwarf.gaia_data[table]['parallax']

        bins = np.linspace(parallax.min(), parallax.max(), num=30)
        ax.hist(parallax, bins=bins)

        ax.set_title(r'radius=${table}^\circ$')
        ax.set_xlabel('Parallax [mas]')
        ax.set_ylabel('Bin counts')

    fig.suptitle('Parallax histograms for different radii')
    fig.savefig(f'{dwarf.path}/plots/parallax_plot.pdf', bbox_inches='tight')


def quiver_plot(dwarf):
    """Create quiver plot of stellar position and proper motion."""
    fig, axs = plt.subplots(ncols=2, nrows=len(dwarf.gaia_data)//2 + len(dwarf.gaia_data) % 2)

    for ax, table in zip(axs, dwarf.gaia_data):
        ra = dwarf.gaia_data[table]['ra']
        dec = dwarf.gaia_data[table]['dec']
        pmra = dwarf.gaia_data[table]['pmra']
        pmdec = dwarf.gaia_data[table]['pmdec']

        pm_mag = np.hypot(pmra, pmdec)
        pmra = pmra / pm_mag
        pmdec = pmdec / pm_mag

        cmap = cm.viridis
        arrows = ax.quiver(ra, dec, pmra, pmdec, color=cmap(pm_mag), clim=(0, pm_mag.max()), units='xy', pivot='tail', width=0.01, headwidth=2, headlength=3, minlength=0.01)
        cbar = colorbar_for_subplot(fig, ax, cmap, image=arrows)
        ax.set_facecolor('xkcd:black')

        ax.set_title(f'radius={table}, density={np.round(len(ra)/(np.pi * table**2), 2)}')

        ax.set_xlabel(r"Right ascension [$^\circ$]")
        ax.set_ylabel(r"Declination [$^\circ$]")
        cbar.ax.set_ylabel("Proper motion magnitude [mas/yr]", rotation=270, labelpad=10)

    fig.suptitle('Quiver plots for different radii')
    fig.savefig(f'{dwarf.path}/plots/quiver_plot.pdf', bbox_inches='tight')


def pm_histogram(dwarf):
    """Create 2d histogram of proper motions from GAIA search."""
    bound = 5
    bins = np.linspace(-bound, bound, num=20*bound)

    fig, axs = plt.subplots(ncols=2, nrows=len(dwarf.gaia_data)//2 + len(dwarf.gaia_data) % 2, sharex=True, sharey=True)

    for ax, table in zip(axs, dwarf.gaia_data):
        pmra = dwarf.gaia_data[table]['pmra']
        pmdec = dwarf.gaia_data[table]['pmdec']

        counts, _, im = ax.hist2d(pmra, pmdec, bins=(bins, bins), vmin=0, cmap='gnuplot')
        title = f'radius={table}, max count={max(counts)}'

        ax.set_title(title)
        ax.set_xlabel(r"Right ascension proper motion [mas/yr])")
        ax.set_ylabel(r"Declination proper motion [mas/yr]")

        cbar = colorbar_for_subplot(fig, ax, cm.gnuplot, image=im)
        cbar.ax.set_ylabel("Bin counts", rotation=270, labelpad=10)

    fig.suptitle('Proper motion histogram for different radii')
    fig.savefig(f'{dwarf.path}/plots/pmhisto_plot.pdf', bbox_inches='tight')


if __name__ == '__main__':
    gaia_search(90, 90, 'test', dump_to_file=False)
    for i in range(16):
        print(random_cones_outside_galactic_plane())
