#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot dwarf data as histograms.

Author: Jack Runburg
Date: 12-06-2019 12:12

"""

import numpy as np
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm


def trim_axes(axes, N):
    """Massage the axes list to have correct length..."""
    axes = axes.flat
    for ax in axes[N:]:
        ax.remove()
    return axes[:N]


def weighted_average(listofvalues):
    """Compute weighted average of list (a,b) with 'a' as the weights and 'b' as the values."""
    sum = 0
    weights = 0
    for [w, v] in listofvalues:
        sum += w*v
        weights += w
    return sum/weights


def weighted_median(values):
    """Compute the weighted median and return the +/- sigma values."""
    low_jvals, ave_jvals, up_jvals = [], [], []
    temp_j = np.asarray(values)
    weightsum = temp_j[:, 0].sum()
    temp_j = temp_j[temp_j[:, 1].argsort()]
    temp_j[:, 0] = np.cumsum(temp_j, axis=0)[:, 0]
    temp_j = np.divide(temp_j, np.array([weightsum, 1]))
    low_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5-.34, side='left')][1]
    ave_jvals = temp_j[np.searchsorted(temp_j[:, 0], 0.5, side='left')][1]
    up_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5+.34, side='right')][1]

    return np.array([low_jvals, ave_jvals, up_jvals])


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


dwarf_specs = ascii.read('./gaia_data/dwarf_info.ecsv', format='ecsv')

with np.load('dwarf_vels.npz', 'rb', allow_pickle=True) as infile:
    # vels = infile['vels']
    pmra = infile['pmra']
    pmdec = infile['pmdec']

plt.close('all')
mpl.rcParams['axes.labelsize'] = 'x-large'

figsize = (10, 8)
bin_width = 2
d = len(dwarf_specs)
rows = 3
cols = d//rows + 1
fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)
axs = trim_axes(axs, d)

for i, ax in enumerate(fig.axes):
    # pmra[i], pmdec[i] = reject_outliers(pmra[i], pmdec[i])
    ramax = pmra[i].max()
    ramin = pmra[i].min()
    decmax = pmdec[i].max()
    decmin = pmdec[i].min()
    # rabins = np.logspace(-2, 1, num=30)
    # decbins = np.logspace(-2, 1, num=30)
    rabins = np.logspace(-2, 1, num=30)
    rabins = np.concatenate((np.flip(np.negative(rabins)), rabins))
    decbins = np.logspace(-2, 1, num=30)
    decbins = np.concatenate((np.flip(np.negative(decbins)), decbins))
    # rabins = np.logspace(-2, 1, num=int((ramax-ramin)/bin_width))
    # decbins = np.logspace(-2, 1, num=int((decmax-decmin)/bin_width))
    # rabins = np.linspace(ramin, ramax, num=int((ramax-ramin)/bin_width))
    # decbins = np.linspace(decmin, decmax, num=int((decmax-decmin)/bin_width))
    # histo, binx, biny = np.histogram2d(pmra[i], pmdec[i], bins=(rabins, decbins))
    counts, xedges, yedges, im = ax.hist2d(pmra[i], pmdec[i], bins=(rabins, decbins), vmin=0, vmax=8)
    ax.set_title(dwarf_specs['MAIN_ID'][i].strip('NAME'))
    ax.set_xscale('symlog', linthreshx=0.1)
    ax.set_yscale('symlog', linthreshy=0.1)
    # plt.colorbar(im, ax=ax)
    # ax.text(bins_edges[len(bins_edges)//4], int(.75 * max(hist)), 'simbad_pm=' + str(round(dwarf_specs['PM_MAG'][i-1], 3)) + '\nmax=' + str(round(bins_edges[np.argmax(hist)], 3)))

# big axes
fig.add_subplot(111, frameon=False)
plt.tight_layout()
fig.colorbar(im, ax=axs.ravel().tolist())
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("Proper motion, right ascension [mas/yr]")
plt.ylabel("Proper motion, declination [mas/yr]")

fig.savefig('./plots/known_dwarf_histos.pdf', bbox_inches='tight')


# random plot
dwarf_specs = ascii.read('./gaia_data/randomcone_info.ecsv', format='ecsv')

with np.load('randomcone_vels.npz', 'rb', allow_pickle=True) as infile:
    # vels = infile['vels']
    pmra = infile['pmra']
    pmdec = infile['pmdec']

figsize = (10, 8)
bin_width = 2
d = len(dwarf_specs)
cols = 3
rows = d//cols + 1
fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=figsize)
axs = trim_axes(axs, d)

for i, ax in enumerate(fig.axes):
    # pmra[i], pmdec[i] = reject_outliers(pmra[i], pmdec[i])
    ramax = pmra[i].max()
    ramin = pmra[i].min()
    decmax = pmdec[i].max()
    decmin = pmdec[i].min()
    rabins = np.logspace(-2, 1, num=30)
    rabins = np.concatenate((np.flip(np.negative(rabins)), rabins))
    decbins = np.logspace(-2, 1, num=30)
    decbins = np.concatenate((np.flip(np.negative(decbins)), decbins))
    # rabins = np.logspace(-2, 1, num=int((ramax-ramin)/bin_width))
    # decbins = np.logspace(-2, 1, num=int((decmax-decmin)/bin_width))
    # rabins = np.linspace(ramin, ramax, num=int((ramax-ramin)/bin_width))
    # decbins = np.linspace(decmin, decmax, num=int((decmax-decmin)/bin_width))
    # histo, binx, biny = np.histogram2d(pmra[i], pmdec[i], bins=(rabins, decbins))
    counts, xedges, yedges, im = ax.hist2d(pmra[i], pmdec[i], bins=(rabins, decbins), vmin=0, vmax=8)
    # ax.set_title('RA: ' + str(np.round(dwarf_specs['RA'][i], 2)) + ', DEC: ' + str(np.round(dwarf_specs['DEC'][i], 2)) + ', RAD: ' + str(np.round(dwarf_specs['GALDIM_MAJAXIS'][i]/2, 2)))
    ax.set_title('RA: ' + str(np.round(dwarf_specs['RA'][i], 2)) + ', DEC: ' + str(np.round(dwarf_specs['DEC'][i], 2)))
    ax.set_xscale('symlog', linthreshx=0.1)
    ax.set_yscale('symlog', linthreshy=0.1)
    # plt.colorbar(im, ax=ax)

# big axes
fig.add_subplot(111, frameon=False)
plt.tight_layout()
fig.colorbar(im, ax=axs.ravel().tolist())
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("Proper motion, right ascension [mas/yr]")
plt.ylabel("Proper motion, declination [mas/yr]")
fig.savefig('./plots/random_cone_histos.pdf', bbox_inches='tight')
plt.show()
