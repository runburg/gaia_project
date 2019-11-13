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
import glob


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


def all_sky():
    """Create all sky plot of dwarf candidates."""
    from mw_plot import MWSkyMap
    from astropy import units as u

    # setup a MWSkyMap instance
    fig = MWSkyMap(projection='hammer')
    # fig, axs = plt.subplots()

    # alpha value for the milkyway image
    fig.imalpha = 1.

    # setup colormap
    fig.cmap = 'jet'

    colors = ['xkcd:mauve', 'xkcd:coral', 'xkcd:pinkish purple', 'xkcd:tangerine', 'xkcd:vermillion', 'xkcd:tomato']

    known = np.loadtxt('/Users/runburg/github/gaia_project/dsph_search/the_search/tuning/tuning_known_dwarfs.txt', delimiter=',', dtype=str)

    ra_known = known[:, 1].astype(np.float)
    dec_known = known[:, 2].astype(np.float)

    fig.s = 20
    fig.mw_scatter(ra_known*u.degree, dec_known*u.degree, 'xkcd:light grey blue')

    for color, file in zip(colors, glob.glob('region_candidates/*.txt')):
        candidate_list = np.loadtxt(file, delimiter=" ")
        file = file.split('_')
        ra = file[1].strip('ra').astype(np.float)/100
        dec = file[2].strip('dec').astype(np.float)/100
        radius = file[3].strip('rad').astype(np.float)/100

        circle_points = get_points_of_circle(ra, dec, radius)

        # use mw_scatter instead of scatter because we want a background
        fig.s = 1
        fig.mw_scatter(circle_points[:, 0]*u.degree, circle_points[:, 1]*u.degree, lighten_color(color, amount=1.5))

        fig.s = 10
        fig.mw_scatter(candidate_list[:, 0]*u.degree, candidate_list[:, 1]*u.degree, color)

    fig.savefig('/Users/runburg/github/gaia_project/dsph_search/all_sky_candidates.pdf')


def get_points_of_circle(ra_center, dec_center, radius):
    """Get coordinates of circle for plotting."""
    n = 200
    ra = np.linspace(0, 360, num=n)
    dec = 90 - np.ones(n)*radius

    coords = [np.array([raa, decc]) for (raa, decc) in zip(ra, dec)]
    coords_prime = np.ones(len(coords))

    for i, (phi, theta) in enumerate(coords):
        # https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        a = spherical_to_cartesian(phi, theta)
        b = spherical_to_cartesian(ra, dec)
        k = np.cross(a, b) / np.linalg.norm(np.cross(a, b))
        rotation_angle = np.arccos(a.dot(b)/np.linalg.norm(a)/np.linalg.norm(b))
        rotated_vector = a + np.sin(rotation_angle)*np.cross(k, a) + (1-np.cos(rotation_angle))*np.cross(k, np.cross(k, a))
        coords_prime[i] = cartesian_to_spherical(rotated_vector)

    return coords_prime


def icrs_to_galactic(ra_icrs, dec_icrs):
    """Return galactic coordinates."""
    from astropy.coordinates import SkyCoord

    coords = SkyCoord(ra_icrs, dec_icrs, unit='deg', frame='icrs')
    return np.array([coords.galactic.b.value, coords.galactic.l.value]).T


def spherical_to_cartesian(ra, dec):
    """Get cartesian values from spherical."""
    ra_rad = ra * np.pi/180
    dec_rad = np.pi/2 - dec * np.pi/180
    return np.array([np.sin(dec_rad)*np.cos(ra_rad), np.sin(dec_rad)*np.sin(ra_rad), np.cos(dec_rad)])


def cartesian_to_spherical(vec):
    """Get spherical values from cartesian."""
    ra_rad = np.arctan(vec[1]/vec[0])
    dec_rad = np.arctan(vec[2]/np.linalg.norm(vec))
    return np.array([ra_rad * np.pi/180, 90 - dec_rad * np.pi/180]).T


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.

    Input can be matplotlib color string, hex string, or RGB tuple. Make amount > 1 to darken.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)

    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except KeyError:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


if __name__ == '__main__':
    all_sky()
