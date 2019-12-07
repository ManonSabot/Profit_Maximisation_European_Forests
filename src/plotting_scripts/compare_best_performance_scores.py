#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to plot the performance scores

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "plot the performance scores of the various configurations"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (10.02.2018)"
__email__ = "m.e.b.sabot@gmail.com"

#=======================================================================

import warnings  # ignore these warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# import general modules
import os  # check for files, paths
import sys  # check for files, paths
import numpy as np  # array manipulations, math operators
import numpy.ma as ma  # masked arrays
from scipy.stats import percentileofscore as psc  # quantile ranks

# plotting modules
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
import string  # automate subplot lettering

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.Utils import get_main_dir  # locate data
from TractLSM.Utils import read_csv  # read in data files

# local modules
from plot_utils import get_best_kmax_calib, get_best_kmax_climate


#=======================================================================

def main(project, sites, split_years):

    """
    Main: calculates the yearly performance scores between April - Nov
          and plots them.

    Arguments:
    ----------
    project: string
        name of the repository which contains the performance file

    sites: array
        names of the sites to plot

    split_years: array
        how to group years, e.g. drought years are 2003 and 2006

    Returns:
    --------
    'performance_scores.png' in output/figures/final_4_paper/

    """

    # paths
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    basedir = os.path.join(basedir, 'output')
    figdir = os.path.join(os.path.join(basedir, 'figures'), 'final_4_paper')

    basedir = os.path.join(os.path.join(basedir, 'projects'), project)
    perfile = os.path.join(basedir, 'perf_scores.csv')

    if not os.path.isdir(figdir):  # create figdir if it doesn't exist
        os.makedirs(figdir)

    # predefining this info for convenience
    configs = ['_', 'H_', 'L_']
    labels = ['Climate', 'Control']

    # look for which is the "best" kmax-calib
    kmax = get_best_kmax_calib(sites, basedir)

    # look for which is the "best" kmax-climate
    climate = get_best_kmax_climate(sites, basedir)

    plot_ranks(perfile, figdir, sites, split_years, labels, configs, kmax=kmax,
               climate=climate)

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def get_scores(fname, configs, sites=None, kmax=None, climate=None,
               split_years=None):

    """
    Retrieves each of the yearly metrics and quantile ranks for all the
    sites and fluxes.

    Arguments:
    ----------
    fname: string
        name (path) of the performance file

    configs: array
        kmax optimisation configurations: high, optimal, low

    sites: array
        site names

    kmax: array
        kmax numbers for the calibration

    climate: array
        which is the best climate for each site

    split_years: array
        how to group years, e.g. drought years are 2003 and 2006

    Returns:
    --------
    NMSE: array
        root mean square errors for each of the variables

    MAE: array
        mean absolute error

    SD: array
        standard deviation

    P5: array
        5th percentile

    P95: array
        95th percentile

    ranks: array
        quantile ranks

    sites: array
        site names

    """

    # read data
    df, __ = read_csv(fname)
    df = df[df['year'] != 'overall']

    # site list
    if sites is None:
        sites = sorted(list(np.unique(df['Site'])))

    if kmax is not None:
        if climate is not None:
            types = ['sample', 'adjust', 'control']

        else:
            types = ['sample', 'adjust_average', 'adjust_extreme', 'control']

    else:
        if climate is not None:
            types = ['adjust', 'control']

        else:
            types = ['adjust_average', 'adjust_extreme', 'control']

    # empty arrays to store the performance scores in
    # dims: sites, types, configs, An or ET, drought or not, years
    NMSE = np.zeros((len(sites), len(types), len(configs), 2, 2, 2))
    MAE = np.zeros((len(sites), len(types), len(configs), 2, 2, 2))
    SD = np.zeros((len(sites), len(types), len(configs), 2, 2, 2))
    P5 = np.zeros((len(sites), len(types), len(configs), 2, 2, 2))
    P95 = np.zeros((len(sites), len(types), len(configs), 2, 2, 2))

    # loop over ET and An
    for i in range(2):

        # loop over configurations (i.e. opt, high, low)
        for j in range(len(configs)):

            # loop over sites
            for k in range(len(sites)):

                for l in range(len(types)):

                    subset = df[df['fname'].str.contains('%s%s' %
                                                         (sites[k],
                                                          configs[j]))]

                    if 'sample' in types[l]:
                        if kmax is not None:
                            subset = df[df['fname'].str.contains('%s%s_' %
                                        (sites[k], kmax[k]))]

                    if 'adjust' in types[l]:
                        if climate is not None:
                            check = (subset['project']
                                     .str.contains(climate[k].lower()))

                            if any(check):
                                subset = subset[check]

                    sub = subset[subset['project'].str.contains(types[l])]

                    for m in range(len(sub)):

                        ssub = sub.iloc[m]

                        if split_years is not None:
                            if any(yr in ssub['fname'] for yr in split_years):
                                n = 0  # drought or not?

                            else:
                                n = 1

                        else:
                            n = 0

                        if i == 0:  # An
                            if np.isclose(NMSE[k, l, j, i, n, 0], 0.):
                                o = 0  # years

                            else:
                                o = 1

                            NMSE[k, l, j, i, n, o] = ssub['NMSE(An)']
                            MAE[k, l, j, i, n, o] = ssub['MAE(An)']
                            SD[k, l, j, i, n, o] = ssub['SD(An)']
                            P5[k, l, j, i, n, o] = ssub['P5(An)']
                            P95[k, l, j, i, n, o] = ssub['P95(An)']

                        else:  # ET
                            if np.isclose(NMSE[k, l, j, i, n, 0], 0.):
                                o = 0

                            else:
                                o = 1

                            NMSE[k, l, j, i, n, o] = ssub['NMSE(ET)']
                            MAE[k, l, j, i, n, o] = ssub['MAE(ET)']
                            SD[k, l, j, i, n, o] = ssub['SD(ET)']
                            P5[k, l, j, i, n, o] = ssub['P5(ET)']
                            P95[k, l, j, i, n, o] = ssub['P95(ET)']

    # deal with configs that don't apply (e.g. 'high' & 'low' to sample)
    for metric in [NMSE, MAE, SD, P5, P95]:
        metric[:, types.index('sample'), 1:, :, :, :] = 0.
        metric[:, types.index('control'), 1:, :, :, :] = 0.

    # mask zeros in array (e.g. only one year, calib, control)
    NMSE = ma.masked_where(NMSE == 0., NMSE)
    MAE = ma.masked_where(MAE == 0., MAE)
    SD = ma.masked_where(SD == 0., SD)
    P5 = ma.masked_where(P5 == 0., P5)
    P95 = ma.masked_where(P95 == 0., P95)

    # empty array to store the ranks in
    # dims: sites, types, configs, An or ET, drought or not, metrics
    ranks = np.zeros((len(sites), len(types), len(configs), 2, 2, 5))

    # loop over ET and An
    for i in range(2):

        # loop over sites
        for k in range(len(sites)):

            # loop over event types, i.e. drought or not
            for n in range(2):

                # loop over the metrics of performance
                for m, metric in enumerate([NMSE, MAE, SD, P5, P95]):

                    # average the metric within an event category
                    a = ma.mean(metric[k, :, :, i, n, :], axis=2)
                    b = a[~ma.getmask(ma.mean(metric[k, :, :, i, n, :],
                                              axis=2))]

                    # ranks computed in weak quantile regression
                    rk = np.array([psc(b, c, 'weak') / 100. for c in b])

                    track = 0  # tracker for the computed ranks

                    # loop over experiment types
                    for l in range(len(types)):

                        # loop over configurations
                        for j in range(len(configs)):

                            # ranks where the metrics are not nan!
                            if a[l, j] == b[track]:
                                ranks[k, l, j, i, n, m] = rk[track]
                                track += 1  # update the tracker

                                if track == len(b):  # all is filled
                                    break

    # average ranking across metrics
    ranks = ma.masked_where(ranks <= 0., ranks)
    ranks = ma.mean(ranks, axis=5)

    return NMSE, MAE, SD, P5, P95, ranks, sites


def plot_ranks(fname, figdir, sites, dryrs, labels, configs, kmax=None,
               climate=None):

    """
    Plots the quantile ranks for all the sites and fluxes, as well as
    summaries across sites.

    Arguments:
    ----------
    fname: string
        name (path) of the performance file

    figdir: string
        name of the repository in which to output the figure

    sites: array
        site names

    dryrs: array
        how to group years, e.g. drought years are 2003 and 2006

    labels: array
        labels for the configuration shown in the plot

    configs: array
        kmax optimisation configurations: high, optimal, low

    kmax: array
        kmax numbers for the calibration

    climate: array
        which is the best climate for each site

    Returns:
    --------
    'performance_scores.png' in output/figures/final_4_paper/

    """

    # plotting colors and marker shape
    colours = ['#016c59', '#762a83']

    if kmax is not None:
        labels = ['Calibration'] + labels
        colours = ['#7fbc41'] + colours

    if climate is not None:  # which climate colours
        ccolours = ['#1f78b4' if e == 'Average' else '#fc4e2a' for e in
                    climate]

    mark = 'o'
    msize = 2.8

    # site indexes
    site_idxs = np.linspace(0., len(sites), len(sites))

    # aspect ratio
    plt.figure(figsize=(4., 2.))
    rspan, cspan = 1, 8

    # declare axes
    ax1 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (0, 5), rowspan=rspan,
                           colspan=cspan)
    ax2 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (rspan, 5),
                           rowspan=rspan, colspan=cspan)
    ax3 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (0, cspan + 5),
                           rowspan=rspan, colspan=cspan)
    ax4 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (rspan, cspan + 5),
                           rowspan=rspan, colspan=cspan)
    ax5 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (0, 0),
                           rowspan=rspan, colspan=3)
    ax6 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (rspan, 0),
                           rowspan=rspan, colspan=3)
    ax7 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (0, cspan * 2 + 6),
                           rowspan=rspan, colspan=3)
    ax8 = plt.subplot2grid((rspan * 2, 2 * cspan + 7), (rspan, cspan * 2 + 6),
                           rowspan=rspan, colspan=3)

    # share the relevant axes
    ax2.get_shared_y_axes().join(ax2, ax4)
    ax1.get_shared_y_axes().join(ax1, ax5)
    ax2.get_shared_y_axes().join(ax2, ax6)
    ax5.get_shared_x_axes().join(ax5, ax6)

    # retrieve the ranks
    __, __, __, __, __, ranks, __ = get_scores(fname, configs, sites=sites,
                                               kmax=kmax, climate=climate,
                                               split_years=dryrs)

    # plot all
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:

        # share the relevant x axes
        if (ax != ax1) and (ax != ax5) and (ax != ax6):
            ax1.get_shared_x_axes().join(ax1, ax)

        boxes = []  # summary box plots

        for i in range(len(labels)):  # go over labels

            if (ax == ax1) or (ax == ax3) or (ax == ax5):
                k = 0

            else:
                k = 1

            if (ax == ax1) or (ax == ax2):
                l = 0

            else:
                l = 1

            if kmax is None:
                ii = i

            else:
                ii = i - 1

            isite = site_idxs + ii / 4.

            # plot the optimal and control ranks
            if (ax == ax1) or (ax == ax2) or (ax == ax3) or (ax == ax4):
                if labels[i] == 'Climate':
                    ax.scatter(isite, ranks[:, i, 0, k, l], color=ccolours,
                               marker=mark, s=msize, label=labels[i])

                else:
                    ax.scatter(isite, ranks[:, i, 0, k, l], color=colours[i],
                               marker=mark, s=msize, label=labels[i])

            elif (ax == ax5) or (ax == ax6):
                spread = (ranks[:, i, 0, k, :]
                               [~ma.getmask(ranks[:, i, 0, k, :])])
                spread = [e for e in spread.flatten() if e]
                boxes += [spread]

                if ((kmax is None) and (i == 1)) or ((kmax is not None) and
                                                     (i == 2)):
                    bp = ax.boxplot(boxes, medianprops=dict(linewidth=0.),
                                    meanprops=dict(linestyle='-',
                                                   linewidth=0.75,
                                                   color='#ffffbf'),
                                    flierprops=dict(marker='o',
                                                    markersize=0.3),
                                    bootstrap=1000, whis=[10, 90], widths=0.4,
                                    patch_artist=True, meanline=True,
                                    showmeans=True)

                    track = 0

                    for box in bp['boxes']:

                        box.set(color='k', linewidth=0.2)
                        box.set(facecolor=colours[track])
                        track += 1

                    track = 0

                    for flyer in bp['fliers']:

                        flyer.set(markeredgecolor=colours[track])
                        track += 1

                    track = 0

                    for whisker in bp['whiskers']:

                        whisker.set(color=colours[int(track / 2)],
                                    linewidth=0.5)
                        track += 1

                    track = 0

                    for cap in bp['caps']:

                        cap.set(color=colours[int(track / 2)],
                                linewidth=0.5)
                        track += 1

            for j in range(len(sites)):  # plot range in behaviours

                if labels[i] == 'Climate':
                    if ((ax == ax1) or (ax == ax3) or (ax == ax2) or
                       (ax == ax4)):
                        ax.vlines(isite[j], ma.amin(ranks[j, i, :, k, l]),
                                  ma.amax(ranks[j, i, :, k, l]),
                                  colors=ccolours[j],
                                  linewidth=0.5, zorder=-1)

    # ax labels and format
    start, end = ax1.get_xlim()
    track = 0

    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]:

        ax.yaxis.set_minor_locator(plt.NullLocator())

        # draw site separation
        if (ax == ax1) or (ax == ax2) or (ax == ax3) or (ax == ax4):

            for isite in site_idxs[:-1]:

                ax.axvline(x=isite + 5. / 9., color='k', alpha=0.25,
                           linestyle=':', linewidth=0.5)

        if ax == ax5:
            left, right = ax.get_xlim()
            left -= 0.25
            right += 0.25

        # deal with y-axis labels
        if (ax == ax3) or (ax == ax4) or (ax == ax5) or (ax == ax6):
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.tick_params(axis='y', length=0)
            ax.set_yticklabels([])

        # format y-tick labels
        if (ax == ax1) or (ax == ax2):
            ax.yaxis.set_major_locator(plt.MaxNLocator(4))
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            for tick in ax.yaxis.get_major_ticks():

                tick.label.set_fontsize(6)

        # deal with x-axis labels and limits
        ax.set_xticklabels([])
        ax.get_xaxis().set_ticks([])
        ax.set_xlim(start, end)

        if (ax == ax5) or (ax == ax6) or (ax == ax7) or (ax == ax8):
            ax.set_xlim(left, right)
            ax.axis('off')

        # arrange site names
        if (ax == ax2) or (ax == ax4):
            ax.set_xticks(site_idxs + 0.1)
            ax.tick_params(axis='x', color='w', length=0.15)
            ax.set_xticklabels(sites, rotation=90, ha='center', fontsize=6.5)

        # add skill arrow
        if (ax == ax7) or (ax == ax8):
            ax.text(1., 0.875, 'low\nskill', va='center', ha='center',
                    fontsize=5.75)
            ax.annotate('high\nskill', xy=(1., 0.75), xytext=(1., 0.1),
                        arrowprops=dict(arrowstyle='<-', linewidth=0.75),
                        va='center', ha='center', fontsize=5.75)

        # indicate years
        if ax == ax1:
            ax.set_title('2003 and 2006', loc='center', fontsize=7)

        if ax == ax3:
            ax.set_title('2002 and 2005', loc='center', fontsize=7)

        if ax == ax5:
            ax.set_title('Overall', loc='center', fontsize=7)

        if track < 4:  # subplot numbering
            ax.text(0.89, 0.89, '(%s)' % (string.ascii_lowercase[track]),
                    transform=ax.transAxes, fontsize=6, fontweight='heavy')

        track += 1

    ax1.set_ylabel('GPP q-rank', labelpad=40, fontsize=7)
    ax2.set_ylabel('ET q-rank', labelpad=40, fontsize=7)

    # adjust figure size
    plt.subplots_adjust(wspace=0.15, hspace=0.075)

    # legend
    if climate is not None:  # plot a text box
        extra_labels = np.unique(climate)
        ccolours = ['#1f78b4', '#fc4e2a']

        for i, label in enumerate(labels):

            pos = 0.02 - i / 6.5

            if label == 'Climate':
                pos2 = pos + 0.0525

                for j, lab in enumerate(extra_labels):

                    ax6.text(3., pos2, lab, {'color': ccolours[j],
                                             'fontsize': 5.75})

                    if j == 0:
                        pos2 -= 0.105

                ax6.text(2.55, pos, r'$\{$', {'color': colours[i],
                                              'fontsize': 7.})

            ax6.text(0., pos, label, {'color': colours[i], 'fontsize': 5.75})

    else:  # plot a legend
        leg = ax6.legend(bbox_to_anchor=(-0.375, -0.0015), loc=2,
                         handlelength=0, handletextpad=0, labelspacing=0.3,
                         fancybox=True, fontsize=5.75, ncol=1)

        track = 0

        for text in leg.get_texts():

            plt.setp(text, color=colours[track], fontweight='bold')
            track += 1

        for handle in leg.legendHandles:

            handle.set_visible(False)

    # save the figure
    namefig = os.path.join(figdir, 'performance_scores')
    plt.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight')
    plt.savefig('%s.eps' % (namefig), dpi=600, bbox_inches='tight')
    plt.close()

    return


#=======================================================================

if __name__ == "__main__":

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.weight'] = 'light'
    plt.rcParams['font.size'] = 7

    # user input
    project = 'var_kmax'  # where the perf scores are stored
    sites = ['Hyytiala', 'Soroe', 'Loobos', 'Hesse', 'Parco', 'Puechabon',
             'Rocca1', 'Rocca2', 'ElSaler1', 'Espirra']  # sites to plot
    dryrs = ['2003', '2006']  # years of drought

    try:
        plt.rc('text', usetex=True)

        try:
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
                                                   r'\usepackage{cmbright}']
            main(project, sites, dryrs)

        except Exception as e:
            try:
                plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
                main(project, sites, dryrs)

            except Exception as e:
                plt.rcParams['text.latex.preamble'] = \
                    [r'\usepackage{cmbright}']
                main(project, sites, dryrs)

    except Exception as e:
        main(project, sites, dryrs)
