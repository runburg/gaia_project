#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make plots for dsph search project.

Author: Jack Runburg
Date: 11-09-2019 11:16


"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl


def colorbar_for_subplot(fig, axs, cmap, image):
    """Place a colorbar by each plot."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(axs)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(image, cax=cax)
    # cm.ScalarMappable(norm=None, cmap=cmap),

    return cbar


def trim_axes(axes, N):
    """Trim the axes list to proper length."""
    if N > 1:
        axes = axes.flat
        for ax in axes[N:]:
            ax.remove()
        return axes[:N]

    return [axes]


def plot_setup(rows, cols, d=0, buffer=(0.3, 0.3)):
    """Set mpl parameters for beautification."""
    # setup plot
    plt.close('all')
    mpl.rcParams['axes.labelsize'] = 'large'
    mpl.rcParams['ytick.labelsize'] = 'x-small'
    mpl.rcParams['xtick.labelsize'] = 'x-small'
    # mpl.rcParams['figure.subplot.wspace'] = buffer[0]
    # mpl.rcParams['figure.subplot.hspace'] = buffer[1]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    figsize = (6 * cols + buffer[0], 5.5 * rows + buffer[1])
    fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)
    if d != 0:
        axs = trim_axes(axs, d)

    return fig, axs


def string_name(dwarf):
    """Give a string rep of name."""
    try:
        ra, _, dec = dwarf.name.partition('_')
        return f'({round(float(ra)/100, 2)}, {round(float(dec)/100, 2)})'
    except ValueError:
        return f'{dwarf.name}'


def parallax_histogram(dwarf):
    """Create a histogram of parallax values."""
    cols = 2
    rows = len(dwarf.gaia_data)//2 + len(dwarf.gaia_data) % 2

    fig, axs = plot_setup(rows, cols, d=len(dwarf.gaia_data))

    for ax, table in zip(axs.flatten(), dwarf.gaia_data):
        parallax = dwarf.gaia_data[table]['parallax']

        bins = np.linspace(parallax.min(), parallax.max(), num=30)
        ax.hist(parallax, bins=bins)

        ax.set_title(rf'radius=${table}^\circ$')
        ax.set_xlabel('Parallax [mas]')
        ax.set_ylabel('Bin counts')

    fig.suptitle(f'Parallax histograms for {string_name(dwarf)}')
    fig.savefig(f'{dwarf.path}/plots/parallax_plot.pdf', bbox_inches='tight')


def quiver_plot(dwarf):
    """Create quiver plot of stellar position and proper motion."""
    cols = 2
    rows = len(dwarf.gaia_data)//2 + len(dwarf.gaia_data) % 2

    fig, axs = plot_setup(rows, cols, d=len(dwarf.gaia_data))

    for ax, table in zip(axs.flatten(), dwarf.gaia_data):
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

    fig.suptitle(f'Quiver plots for {string_name(dwarf)}')
    fig.savefig(f'{dwarf.path}/plots/quiver_plot.pdf', bbox_inches='tight')


def pm_histogram(dwarf):
    """Create 2d histogram of proper motions from GAIA search."""
    bound = 5
    bins = np.linspace(-bound, bound, num=20*bound)

    cols = 2
    rows = len(dwarf.gaia_data)//2 + len(dwarf.gaia_data) % 2

    fig, axs = plot_setup(rows, cols, d=len(dwarf.gaia_data))

    for ax, table in zip(axs.flatten(), dwarf.gaia_data):
        pmra = dwarf.gaia_data[table]['pmra']
        pmdec = dwarf.gaia_data[table]['pmdec']

        counts, _, _, im = ax.hist2d(pmra, pmdec, bins=(bins, bins), vmin=0, cmap='gnuplot')
        title = f'radius={table}, max count={str(counts.max())}'

        ax.set_title(title)
        ax.set_xlabel(r"Right ascension pm [mas/yr])")
        ax.set_ylabel(r"Declination pm [mas/yr]")

        cbar = colorbar_for_subplot(fig, ax, cm.gnuplot, image=im)
        cbar.ax.set_ylabel("Bin counts", rotation=270, labelpad=10)

    fig.suptitle(f'Proper motion histogram for {string_name(dwarf)}')
    fig.savefig(f'{dwarf.path}/plots/pmhisto_plot.pdf', bbox_inches='tight')


def mag_v_color(dwarf):
    """Create color magnitude diagram."""
    cols = 2
    rows = len(dwarf.gaia_data)//2 + len(dwarf.gaia_data) % 2

    fig, axs = plot_setup(rows, cols, d=len(dwarf.gaia_data))

    for ax, table in zip(axs.flatten(), dwarf.gaia_data):
        color = dwarf.gaia_data[table]['bp_rp']
        mag = dwarf.gaia_data[table]['phot_g_mean_mag']

        lines = ax.scatter(color, mag, c=color, cmap=cm.cool)
        title = f'radius={table}, color-mag diagram'

        ax.set_title(title)
        ax.set_xlabel(r"BpRp [mag]")
        ax.set_ylabel(r"photGMeanMag [mag]")

    fig.suptitle(f'Color magnitude diagrams for {string_name(dwarf)}')
    fig.savefig(f'{dwarf.path}/plots/colormag.pdf', bbox_inches='tight')


def all_sky(candidate_list):
    """Create all sky plot of dwarf candidates."""
    from mw_plot import MWSkyMap
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    # setup a MWSkyMap instance
    fig = MWSkyMap(projection='hammer')

    # alpha value for the milkyway image
    fig.imalpha = 1.

    # setup colormap
    fig.cmap = 'jet'

    # use mw_scatter instead of scatter because we want a colorbar
    fig.mw_scatter(candidate_list[:, 0]*u.degree, candidate_list[:, 1]*u.degree, "xkcd:mauve")

    fig.savefig(file='./dsph_search/all_sky_candidates.pdf')


if __name__ == '__main__':
    coords = np.loadtxt("./dsph_search/candidate_coords.txt", delimiter=" ")
    all_sky(coords)
