#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to plot the fluxes by grouping sites by groups of 5

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "plot the site-level yearly fluxes of the configurations"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (10.02.2018)"
__email__ = "m.e.b.sabot@gmail.com"

#=======================================================================

import warnings  # ignore these warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# import general modules
import argparse  # read in the user input
import os  # check for files, paths
import sys  # check for files, paths
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files

# plotting modules
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator  # tick locators
from matplotlib.ticker import FormatStrFormatter  # tick labels
import string  # automate subplot lettering
from matplotlib.lines import Line2D  # custom legend

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.Utils import get_main_dir  # locate data
from TractLSM.Utils import read_csv  # read in data files

# local plotting modules
from plot_utils import find_all_combis, find_site_data  # locate
from plot_utils import get_best_kmax_calib  # best kmax's number name
from plot_utils import convert_units  # units conversions
from plot_utils import split  # split arrays
from plot_utils import select_dates, running_mean, arrays_boundary
from plot_utils import plot_var_2_axis, add_bars_2_axis  # plot


#=======================================================================

def main(project1, project2, sites, split_years):

    """
    Main: finds the "best" performing configurations and plots them.

    Arguments:
    ----------
    project1: string
        name of the repository which contains the modelled data

    project1: string
        name of the repository which contains the modelled data

    sites: array
        names of the sites to plot

    split_years: array
        how to group years, e.g. drought years are 2003 and 2006

    Returns:
    --------
    files of the form 'what_north/south_drought/no_drought.png' in
    output/figures/final_4_paper/

    """

    pd.options.mode.chained_assignment = None  # disable these warnings

    # paths
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    basedir = os.path.join(basedir, 'output')
    figdir1 = os.path.join(os.path.join(basedir, 'figures'), 'final_4_paper')
    figdir2 = os.path.join(os.path.join(basedir, 'figures'),
                           'not_shown_in_paper')
    basedir = os.path.join(os.path.join(basedir, 'projects'), 'var_kmax')

    if not os.path.isdir(figdir1):  # create figdir if it doesn't exist
        os.makedirs(figdir1)

    if not os.path.isdir(figdir2):  # create figdir if it doesn't exist
        os.makedirs(figdir2)

    if 'sample' in project1:  # look for which is the "best" kmax
        kmax = get_best_kmax_calib(sites, basedir)

    elif project2 is not None:
        if 'sample' in project2:
            kmax = get_best_kmax_calib(sites, basedir)

        else:
            kmax = None

    else:
        kmax = None

    if project2 is None:  # no comparison, only plot the project
        if kmax is None:
            fluxes(figdir2, sites, split_years, project1)

        else:
            fluxes(figdir1, sites, split_years, project1, kmax=kmax)

    else:  # comparison
        fluxes(figdir1, sites, split_years, project1, project2=project2,
               kmax=kmax)  # best configurations
        fluxes(figdir2, sites, split_years, project1, project2=project2,
               kmax=kmax, which='worst')  # worst configurations

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def each_flux(site, year, project, ax1, ax2, fill_between=True,
              right_labels=True):

    """
    Retrieves the data of a specific year at a single site and plots the
    fluxes, LAI, and precipitation.

    Arguments:
    ----------
    site: string
        individual site name

    year: int
        individual year to plot

    project: string
        name of the repository which contains the data to plot

    ax1: matplotlib object
        axis on which the data must be plotted

    ax2: matplotlib object
        axis on which the data must be plotted

    fill_between: boolean
        if True, ranges between the lowest and the highest version of
        the configuration will be plotted. This is useful for the kmax
        optimal behaviours.

    right_labels: boolean
        if True, ticks and labels will be added on the right side for
        the secondary y axes.

    Returns:
    --------
    Plots the relevant data on ax1 and ax2.

    """

    # declare axes
    ax3 = ax1.twinx()  # LAI axis
    ax4 = ax2.twinx()  # precip axis

    # make sure the twin axes are behind the original axes
    ax3.set_zorder(ax1.get_zorder() - 1)
    ax4.set_zorder(ax2.get_zorder() - 1)

    # number of ticks per axis
    ax1.yaxis.set_major_locator(MaxNLocator(4, integer=True, prune='lower'))
    ax2.yaxis.set_major_locator(MaxNLocator(4, integer=True, prune='lower'))
    ax3.yaxis.set_major_locator(MaxNLocator(3, prune='lower'))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # colours (grey, blue, purple, green)
    col = ['dimgrey', '#4393c3', '#762a83', '#7fbc41']

    # running mean
    N = 14  # fortnightly

    # find all related data
    idata, odata = find_site_data(site, year, project)

    # specifics
    if ((any('%sH' % (site) in e for e in idata) and
         any('%sL' % (site) in e for e in idata)) and not
       any('%sk%d' % (site, str.isdigit(e)) in e for e in idata)):
        order = ['H', '', 'L']
        labels = ['high', 'opt.', 'low']

    else:
        fill_between = False  # reset to false
        order = ['', ] * len(idata)
        labels = ['calib.', ] * len(idata)

    # plot the obs
    for data in idata:

        if ('%s_' % (site)) in data:

            # retrieve and convert the data
            df, __ = read_csv(data)
            delta = df['hod'].iloc[1] - df['hod'].iloc[0]

            # restrict to day time and acceptable LAI
            mask = np.logical_and(df['PPFD'] > 50., df['LAI'] > 0.001)
            df[['GPP', 'Qle']] = (df[['GPP', 'Qle']].where(mask)
                                                    .fillna(value=0.))
            df = convert_units(df)

            # dates
            start_date = pd.to_datetime(int(year * 1000 + df['doy'].iloc[0]),
                                        format='%Y%j')

            if df['doy'].iloc[-1] > df['doy'].iloc[0]:
                end_date = pd.to_datetime(int(year * 1000 +
                                              df['doy'].iloc[-1]),
                                          format='%Y%j')

            else:
                end_date = pd.to_datetime(int((year + 1) * 1000 +
                                          df['doy'].iloc[-1]), format='%Y%j')

            dates = pd.date_range(start=start_date, end=end_date, freq='D')

            # northern hemisphere, start in April, go to Nov
            dates_2_plot = select_dates(df, year)

            # running means on iQle and iGPP
            iGPP = running_mean(df['GPP'], N, delta)
            iQle = running_mean(df['Qle'], N, delta)
            iGPP = iGPP.iloc[dates_2_plot]
            iQle = iQle.iloc[dates_2_plot]
            plot_var_2_axis([iGPP, iQle], [ax1, ax2],
                            dates=dates[dates_2_plot], color='k', label='Obs.',
                            share_axes=True)

            # LAI clim with GPP on plot
            LAI = (df['LAI'].groupby(df.index // int(24. / delta) *
                                     int(24. / delta))
                            .mean()
                            .iloc[dates_2_plot])
            LAI.index = dates[dates_2_plot]
            ax3.plot(LAI.index.to_pydatetime(), LAI, color=col[0],
                     linewidth=0.5)
            ax3.tick_params(axis='y')

            # daily precip from top of y axis on ET plot
            precip = (df['precip'].groupby(df.index // int(24. / delta) *
                                           int(24. / delta))
                                  .mean()
                                  .iloc[dates_2_plot])
            add_bars_2_axis([precip], [ax4], dates=dates[dates_2_plot],
                            color=col[1])
            ax4.tick_params(axis='y')

    plot_ctrl = True

    if fill_between:
        store_GPP_psi = [[], []]
        store_Qle_psi = [[], []]

    for i in range(len(order)):

        for data in odata:

            conf = order[i]

            if conf == 'H':
                ls = '--'

            if conf == '':
                ls = '-'

            if conf == 'L':
                ls = ':'

            if ('%s%s_' % (site, conf)) in data:

                # retrieve and convert the data
                df, __ = read_csv(data)
                df.fillna(value=0., inplace=True)

                if plot_ctrl:
                    df1 = convert_units(df, otype='std')
                    oGPP = df1['GPP']
                    oQle = df1['Qle']
                    oGPP = running_mean(oGPP, N, delta)
                    oQle = running_mean(oQle, N, delta)
                    oGPP = oGPP.iloc[dates_2_plot]
                    oQle = oQle.iloc[dates_2_plot]

                    plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                    dates=dates[dates_2_plot], color=col[2],
                                    label='Control', share_axes=True)
                    plot_ctrl = False

                for j, psi in enumerate(['psi1', 'psi2']):

                    try:
                        df1 = convert_units(df, otype=psi)
                        oGPP = df1['GPP']
                        oQle = df1['Qle']
                        oGPP = running_mean(oGPP, N, delta)
                        oQle = running_mean(oQle, N, delta)
                        oGPP = oGPP.iloc[dates_2_plot]
                        oQle = oQle.iloc[dates_2_plot]

                        if fill_between:  # plot the sim(s)
                            store_GPP_psi[j] += [oGPP]
                            store_Qle_psi[j] += [oQle]

                            if i == len(order) - 1:
                                k = order.index('')
                                label = r'Profit$_{\rm max}$'
                                plot_var_2_axis([store_GPP_psi[j][k],
                                                 store_Qle_psi[j][k]],
                                                [ax1, ax2],
                                                dates=dates[dates_2_plot],
                                                color=col[3], label=label,
                                                share_axes=True)  # plot

                                min_GPP = arrays_boundary(store_GPP_psi[j],
                                                          'min')
                                max_GPP = arrays_boundary(store_GPP_psi[j],
                                                          'max')
                                min_Qle = arrays_boundary(store_Qle_psi[j],
                                                          'min')
                                max_Qle = arrays_boundary(store_Qle_psi[j],
                                                          'max')

                                oGPP.index = dates[dates_2_plot]
                                pydates = oGPP.index.to_pydatetime()
                                ax1.fill_between(pydates, min_GPP, max_GPP,
                                                 facecolor=col[3], alpha=0.6)
                                ax2.fill_between(pydates, min_Qle, max_Qle,
                                                 facecolor=col[3], alpha=0.6)

                        else:
                            label = (r'Profit$_{\rm max}$(k$_{\rm %s}$)' %
                                     (labels[i]))
                            plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                            dates=dates[dates_2_plot],
                                            linestyle=ls, color=col[3],
                                            label=label, share_axes=True)

                    except KeyError:
                        pass

    # precip to top of y axis
    ax4.yaxis.set_label_coords(1.125, 0.7)

    # tighter limits
    bottom, top = ax1.get_ylim()
    ax1.set_ylim(0., top)

    bottom, top = ax2.get_ylim()
    ax2.set_ylim(0., top)

    bottom, top = ax3.get_ylim()

    if bottom < 0.5:  # start LAI at 0 if it ever gets low
        ax3.set_ylim(0., top)

    if not right_labels:
        plt.setp(ax3.get_yticklabels(), visible=False)
        ax3.tick_params(axis='y', length=0)
        plt.setp(ax4.get_yticklabels(), visible=False)
        ax4.tick_params(axis='y', length=0)

    return


def fluxes(figdir, sites, dryrs, project1, project2=None, kmax=None,
           which='best'):

    """
    Retrieves the data of a specific year at a single site and plots the
    fluxes, LAI, and precipitation.

    Arguments:
    ----------
    figdir: string
        name (path) of the repository in which to output the plots

    sites: array
        names of the sites to plot

    dryrs: array
        how to group years, e.g. drought years are 2003 and 2006

    project1: string
        name of the repository which contains the modelled data

    project2: string
        name of the repository which contains the modelled data

    kmax: array
        kmax numbers for the calibration

    which: string
        'best' picks the best performing configurations, 'worst' the
         worst ones

    Returns:
    --------
    files of the form 'what_north/south_drought/no_drought.png' in
    figdir.

    """

    # colours (grey, blue, purple, green)
    col = ['dimgrey', '#4393c3', '#762a83', '#7fbc41']

    # break combis into drought years and non drought years
    allsites, years = find_all_combis(project1)
    csites = [[], []]
    cyears = [[], []]
    csites[0] = [e for i, e in enumerate(allsites) if str(years[i]) in dryrs]
    cyears[0] = [e for e in years if str(e) in dryrs]
    csites[1] = [e for i, e in enumerate(allsites) if not str(years[i]) in
                 dryrs]
    cyears[1] = [e for e in years if not str(e) in dryrs]

    if project2 is not None:
        try:
            basedir = get_main_dir()

            while 'src' in basedir:
                basedir = os.path.dirname(get_main_dir())

            basedir = os.path.join(os.path.join(basedir, 'output'), 'projects')
            fbest = open(os.path.join(os.path.join(basedir,
                                                   'var_kmax_best_climate'),
                                      'best.txt'), 'r')
            flines = fbest.readlines()
            projects = [None] * len(sites)

            for l in flines:

                site = l.split(':')[0].strip()
                pro = l.split(':')[1].strip()
                projects[sites.index(site)] = pro

        except IOError:
            projects = np.tile(project1, len(sites))

    else:
        projects = np.tile(project1, len(sites))

    if which == 'worst':  # do the opposite to the best climate
        projects = [project1 if e != project1 else project2 for e in projects]

    sites = split(sites, 2)
    projects = split(list(projects), 2)

    if kmax is not None:
        kmax = split(kmax, 2)

    # loop over drought years, non drought years
    for dr in range(2):

        north = True

        # loop over northermost, southernmost sites
        for s in range(len(sites)):

            # declare figure
            plt.figure(figsize=(4., 4.))
            rspan = 12
            cspan = 3

            # keep track of things
            track = 0
            Nrow = 0
            Ncol = 0
            right_labels = False

            # loop over sites
            for i in range(len(sites[s])):

                site = sites[s][i]
                project = projects[s][i]
                idx = [ii for ii, e in enumerate(csites[dr]) if e == site]
                subyears = sorted([cyears[dr][ii] for ii in idx])

                if kmax is not None:
                    kn = kmax[s][i]

                # declare axes, depends on number of plots to make
                for j in range(csites[dr].count(site)):

                    year = subyears[j]
                    colspan = cspan * 4

                    if j == 0:
                        ax1 = plt.subplot2grid((rspan * 8 + 3,
                                                cspan * 5 * 3 + 10),
                                               (Nrow, Ncol + j * colspan),
                                               rowspan=rspan, colspan=colspan)
                        ax2 = plt.subplot2grid((rspan * 8 + 3,
                                                cspan * 5 * 3 + 10),
                                               (Nrow + rspan,
                                                Ncol + j * colspan),
                                               rowspan=rspan, colspan=colspan,
                                               sharex=ax1)
                        ax1.text(0.02, 0.85, '(%s)' %
                                 (string.ascii_lowercase[track]),
                                 transform=ax1.transAxes, weight='bold')
                        ax2.text(0.02, 0.85, '(%s)' %
                                 (string.ascii_lowercase[track+1]),
                                 transform=ax2.transAxes, weight='bold')

                        # plot the fluxes
                        if kmax is None:
                            each_flux(site, year, project, ax1, ax2,
                                      right_labels=right_labels)

                        else:
                            each_flux(site + kn, year, project, ax1, ax2,
                                      right_labels=right_labels)

                        plt.setp(ax1.get_xticklabels(), visible=False)
                        ax1.tick_params(axis='x', length=0)

                    if j == 1:
                        ax3 = plt.subplot2grid((rspan * 8 + 3,
                                                cspan * 5 * 3 + 10),
                                               (Nrow, Ncol + j * colspan),
                                               rowspan=rspan, colspan=colspan,
                                               sharey=ax1)
                        ax4 = plt.subplot2grid((rspan * 8 + 3,
                                                cspan * 5 * 3 + 10),
                                               (Nrow + rspan,
                                                Ncol + j * colspan),
                                               rowspan=rspan, colspan=colspan,
                                               sharey=ax2, sharex=ax3)
                        ax3.text(0.02, 0.85, '(%s)' %
                                 (string.ascii_lowercase[track]),
                                 transform=ax3.transAxes, weight='bold')
                        ax4.text(0.02, 0.85, '(%s)' %
                                 (string.ascii_lowercase[track+1]),
                                 transform=ax4.transAxes, weight='bold')

                        # plot the fluxes
                        if kmax is None:
                            each_flux(site, year, project, ax3, ax4)

                        else:
                            each_flux(site + kn, year, project, ax3, ax4)

                        plt.setp(ax3.get_xticklabels(), visible=False)
                        ax3.tick_params(axis='x', length=0)

                        plt.setp(ax3.get_yticklabels(), visible=False)
                        ax3.tick_params(axis='y', length=0)

                        plt.setp(ax4.get_yticklabels(), visible=False)
                        ax4.tick_params(axis='y', length=0)

                    track += 2

                if project2 is not None:
                    if 'extreme' in project:
                        project = 'Extreme'

                    if 'average' in project:
                        project = 'Average'

                    if 'sample' in project:
                        project = 'Calibration'

                    ax1.set_title('%s: %s' % (site, project), loc='left')

                else:
                    ax1.set_title(site, loc='left')

                if Ncol == 0:
                    ax1.set_ylabel('GPP', fontsize=5.5)
                    ax2.set_ylabel('ET', fontsize=5.5)
                    ax1.get_yaxis().set_label_coords(-0.3, 0.5)
                    ax2.get_yaxis().set_label_coords(-0.3, 0.5)

                if track >= 12:
                    if track >= 14:
                        ax1.set_ylabel('LAI', fontsize=5.5)
                        ax2.set_ylabel('PPT', fontsize=5.5)
                        ax1.get_yaxis().set_label_coords(1.45, 0.5)
                        ax2.get_yaxis().set_label_coords(1.45, 0.5)

                    else:
                        ax3.set_ylabel('LAI', fontsize=6)
                        ax4.set_ylabel('PPT', fontsize=6)
                        ax3.get_yaxis().set_label_coords(1.45, 0.5)
                        ax4.get_yaxis().set_label_coords(1.45, 0.5)

                # iterate
                Nrow += 2 * rspan + 4

                if Nrow > rspan * 6:
                    Nrow = 0
                    Ncol += colspan * 2 + 8
                    right_labels = True

                if (track < 12 and Ncol <= 20) or (track > 13 and Ncol > 20):
                    plt.setp(ax1.get_xticklabels(), visible=False)
                    ax1.tick_params(axis='x', length=0.2)

                    if track < 16:
                        plt.setp(ax2.get_xticklabels(), visible=False)
                        ax2.tick_params(axis='x', length=0.2)

                    if (track < 12 and Ncol <= 20):
                        plt.setp(ax3.get_xticklabels(), visible=False)
                        ax3.tick_params(axis='x', length=0.2)
                        plt.setp(ax4.get_xticklabels(), visible=False)
                        ax4.tick_params(axis='x', length=0.2)

            # add legend
            handles, labels = ax2.get_legend_handles_labels()
            handles += [Line2D([0, 1], [1, 1], linestyle='-', color=col[0],
                               lw=0.5),
                        Line2D([0, 1], [1, 1], linestyle='-', color=col[1],
                               lw=0.5)]
            labels += ['LAI', 'PPT']

            ax2.legend(handles, labels, bbox_to_anchor=(-0., -0.5), loc=2,
                       fancybox=True, fontsize=5)

            # adjust figure size
            plt.subplots_adjust(wspace=0.7, hspace=1.2)

            if 'extreme' in project1:
                namefig = 'extreme_'

            if 'average' in project1:
                namefig = 'average_'

            if 'sample' in project1:
                namefig = 'calib_'

            if project2 is not None:
                namefig = '%s_' % (which)

            if north:
                namefig += 'north_'

            else:
                namefig += 'south_'

            if dr == 1:
                if project2 is not None:
                    figdir = figdir.replace('final_4_paper',
                                            'not_shown_in_paper')

                namefig += 'no_'

            namefig = os.path.join(figdir, '%sdrought' % namefig)
            plt.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight',
                        transparent=True)
            plt.savefig('%s.pdf' % (namefig), dpi=1800, bbox_inches='tight',
                        transparent=True)

            north = False

    return


#=======================================================================

if __name__ == "__main__":

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.weight'] = 'light'
    plt.rcParams['font.size'] = 5.
    plt.rcParams['axes.titlepad'] = 1.8

    # user input
    sites = ['Hyytiala', 'Soroe', 'Loobos', 'Hesse', 'Parco', 'Puechabon',
             'Rocca1', 'Rocca2', 'ElSaler1', 'Espirra']  # sites to plot
    dryrs = ['2003', '2006']  # years of drought

    # define the argparse settings
    description = "Plot all the fluxes in 2 figures for a given folder"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('project1', type=str, help='project where output is')
    parser.add_argument('-c', '--project2', type=str, default=None,
                        help='project where output is')
    args = parser.parse_args()

    try:
        plt.rc('text', usetex=True)

        try:
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
                                                   r'\usepackage{cmbright}']
            main(args.project1, args.project2, sites, dryrs)

        except Exception as e:
            try:
                plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
                main(args.project1, args.project2, sites, dryrs)

            except Exception as e:
                plt.rcParams['text.latex.preamble'] = \
                    [r'\usepackage{cmbright}']
                main(args.project1, args.project2, sites, dryrs)

    except Exception as e:
        main(args.project1, args.project2, sites, dryrs)
