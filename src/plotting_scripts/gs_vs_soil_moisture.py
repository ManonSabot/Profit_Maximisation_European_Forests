#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to plot the gs-theta or gs-Ps function shapes

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "plot the site-level overall gs-soil water functions"
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
from pygam import LinearGAM  # fit the functional shapes

# plotting modules
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes   # zoom
from matplotlib.ticker import FormatStrFormatter  # tick labels
import string  # automate subplot lettering
import matplotlib.patheffects as PathEffects  # white shadow around text
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
from plot_utils import split  # split arrays
from plot_utils import select_dates, running_mean  # process data
from plot_utils import plot_var_2_var  # plot


#=======================================================================

def main(project1, project2, sites):

    """
    Main: finds the "best" performing configurations and their leaf
          water to soil water relationships.

    Arguments:
    ----------
    project1: string
        name of the repository which contains the modelled data

    project2: string
        name of the repository which contains the modelled data

    sites: array
        names of the sites to plot

    Returns:
    --------
    files of the form 'what_all/each_water_stress.png' in
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
            soil_moisture_curves(figdir2, sites, project1)
            soil_moisture_curves(figdir2, sites, project1, plot='all')

        else:
            soil_moisture_curves(figdir1, sites, project1, kmax=kmax)
            soil_moisture_curves(figdir1, sites, project1, kmax=kmax,
                                 plot='all')

    else:  # comparison
        soil_moisture_curves(figdir2, sites, project1, project2=project2,
                             kmax=kmax)  # best configurations
        soil_moisture_curves(figdir2, sites, project1, project2=project2,
                             kmax=kmax, which='worst')
        soil_moisture_curves(figdir2, sites, project1, project2=project2,
                             kmax=kmax, plot='all')  # best config.
        soil_moisture_curves(figdir2, sites, project1, project2=project2,
                             kmax=kmax, plot='all', which='worst')

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def each_water_stress(site, years, project, ax, colours):

    """
    Retrieves the data for all years at a single site and plots both the
    gs-theta (of predawn) single points and the associated fitted
    functional forms.

    Arguments:
    ----------
    site: string
        individual site name

    years: array
        all years to plot

    project: string
        name of the repository which contains the data to plot

    ax: matplotlib object
        axis on which the data must be plotted

    colours: array
        colours for each of the configurations

    Returns:
    --------
    Plots the relevant data on the axis.

    """

    # find all related sim data
    odata = []
    idata = []

    for year in years:

        iidata, oodata = find_site_data(site, year, project)
        odata += [oodata]
        idata += [iidata]

    # flatten the data arrays
    odata = [e for data in odata for e in data]
    idata = [e for data in idata for e in data]

    # specifics
    if ((any('%sH' % (site) in e for e in idata) and
         any('%sL' % (site) in e for e in idata)) and not
       any('%sk%d' % (site, str.isdigit(e)) in e for e in idata)):
        order = ['H', '', 'L']
        labels = ['high', 'opt.', 'low']

    else:
        order = ['', ]
        labels = ['calib.', ]

    # append all the sims per config in a single df
    dfs = []
    dates_2_plot = []

    for i in range(len(order)):

        append = False

        for data in odata:

            conf = order[i]
            year = int(data.split('_')[-1].split('.csv')[0])

            if ('%s%s_' % (site, conf)) in data:
                if not append:  # first occurence, declare df
                    df, __ = read_csv(data)
                    dates = select_dates(df, year)

                    # expand to all data per day
                    dates = np.repeat(dates, int(24 / (df['hod'].iloc[1] -
                                      df['hod'].iloc[0])))[:-1]
                    append = True

                else:  # append to existing df
                    ndf, __ = read_csv(data)
                    df = pd.concat([df, ndf], ignore_index=True)
                    dates = select_dates(ndf, year)

                    # expand to all data per day
                    dates = np.repeat(dates, int(24 / (df['hod'].iloc[1] -
                                      df['hod'].iloc[0])))[:-1]

                dates_2_plot += [dates]

        dfs += [df]

    # masks for all years
    dates = [e for sublist in dates_2_plot for e in sublist]
    dates = split(dates, len(dfs))

    # plot gs from the sim(s)
    min_x = []
    max_x = []

    for i in range(len(order)):

        df = dfs[i]
        dates_2_plot = dates[i]
        conf = order[i]

        # filter by months
        delta = df['hod'].iloc[1] - df['hod'].iloc[0]
        df.fillna(value=0., inplace=True)
        df = df.iloc[dates_2_plot]
        df.reset_index(inplace=True)
        Dsteps = int(24. / delta)  # how many timesteps in a day?

        # exclude data 48 hours after rain (replace those values by 0.)
        irain = list(df[df['sw(%s)' % ('std')].diff().gt(0) == True].index
                     .values)

        for psi in ['psi1', 'psi2']:

            try:
                irain += list(df[df['sw(%s)' % (psi)].diff().gt(0) == True]
                              .index.values)

            except KeyError:
                pass

        irain = [list(range(idx, idx + 2 * Dsteps)) for idx in
                 np.unique(irain)]
        irain = np.unique([idx for sublist in irain for idx in sublist if
                          (idx in df.index) and (idx < len(df))])
        df.loc[irain, df.columns[df.columns.str.contains('gs')]] = 0.

        if i == 0:
            ogs = df['gs(%s)' % ('std')]
            osw = df['sw(%s)' % ('std')]

            # only keep the gs 3 hours around midday
            ogs = ogs.where(np.logical_and(df['hod'] >= 9., df['hod'] <= 15.),
                            0.)
            ogs = running_mean(ogs, 1, delta, how='mean_non_zero')
            osw = osw.iloc[::Dsteps]  # predawn

            # make sure the lengths match
            ogs = ogs[:len(osw)]
            osw = osw[:len(ogs)]
            osw.reset_index(inplace=True, drop=True)
            ogs.reset_index(inplace=True, drop=True)

            # drop the nans where they appear!
            osw = osw[ogs.notnull()]
            ogs = ogs.dropna()

            # deal with potential outliers
            mask = ogs < ogs.quantile(.9)
            osw = osw[mask]
            ogs = ogs[mask]
            osw.reset_index(inplace=True, drop=True)
            ogs.reset_index(inplace=True, drop=True)

            # normalised gs
            ogs = ogs / np.amax(ogs.values)

            # 9th order gam fit for the gs-sw function
            gam = (LinearGAM(n_splines=12, spline_order=9,
                             constraints='monotonic_inc')
                   .gridsearch(osw.values.reshape(-1, 1), ogs.values))
            asw = np.linspace(np.amin(osw), np.amax(osw), num=500)
            ags = gam.predict(asw)
            y_int1 = gam.confidence_intervals(asw, width=0.95)

            plot_var_2_var([ogs, ags], [osw, asw], ax, color=colours[0],
                           label='Control')  # plot

            # add the interval confidence to the fit
            ax.fill_between(asw, y_int1[:, 0], y_int1[:, 1], alpha=0.2,
                            facecolor=colours[0], zorder=-1)

            min_x += [np.amin(osw)]
            max_x += [np.amax(osw)]

        for psi in ['psi1', 'psi2']:

            try:
                ogs = df['gs(%s)' % (psi)]
                osw = df['sw(%s)' % (psi)]

                # only keep the gs 3 hours around midday
                ogs = ogs.where(np.logical_and(df['hod'] >= 9.,
                                               df['hod'] <= 15.), 0.)
                ogs = running_mean(ogs, 1, delta, how='mean_non_zero')
                osw = osw.iloc[::Dsteps]  # predawn

                # make sure the lengths match
                ogs = ogs[:len(osw)]
                osw = osw[:len(ogs)]
                osw.reset_index(inplace=True, drop=True)
                ogs.reset_index(inplace=True, drop=True)

                # drop the nans where they appear!
                osw = osw[ogs.notnull()]
                ogs = ogs.dropna()

                # deal with potential outliers
                mask = ogs < ogs.quantile(.9)
                osw = osw[mask]
                ogs = ogs[mask]
                osw.reset_index(inplace=True, drop=True)
                ogs.reset_index(inplace=True, drop=True)

                # normalised gs
                ogs = ogs / np.amax(ogs.values)

                # 9th order gam fit for the gs-sw function
                gam = (LinearGAM(n_splines=12, spline_order=9,
                                 constraints='monotonic_inc')
                       .gridsearch(osw.values.reshape(-1, 1), ogs.values))
                asw = np.linspace(np.amin(osw), np.amax(osw), num=500)
                ags = gam.predict(asw)
                y_int1 = gam.confidence_intervals(asw, width=0.95)

                # plot
                label = r'Profit$_{\rm max}$(k$_{\rm %s}$)' % (labels[i])

                if labels[i] == 'opt.':
                    zorder = -2

                else:
                    zorder = -3

                plot_var_2_var([ogs, ags], [osw, asw], ax, color=colours[1+i],
                               label=label, zorder=zorder)

                # add the interval confidence to the fit
                ax.fill_between(asw, y_int1[:, 0], y_int1[:, 1], alpha=0.2,
                                facecolor=colours[1+i], zorder=zorder)

                min_x += [np.amin(osw)]
                max_x += [np.amax(osw)]

            except KeyError:
                pass

    margin = 0.025 * (np.abs(np.amax(max_x)) - np.abs(np.amin(min_x)))
    ax.set_xlim(np.amin(min_x) - margin, np.amax(max_x) + margin)

    return


def all_water_stress(site, years, project, ax, axin, colour, zorder=0):

    """
    Retrieves the data for all years and all sites and plots each site's
    fitted functional forms of gs-Ps (of predawn).

    Arguments:
    ----------
    site: string
        individual site name

    years: array
        all years to plot

    project: string
        name of the repository which contains the data to plot

    ax: matplotlib object
        axis on which the data must be plotted

    axin: matplotlib object
        inset axis on which zoomed data must be plotted

    colour: array
        colour for the individual site

    zorder: int
        location of the individual site data relative to the other sites

    Returns:
    --------
    Plots the relevant data on the ax.

    """

    # find all related sim data
    odata = []
    idata = []

    for year in years:

        iidata, oodata = find_site_data(site, year, project)
        odata += [oodata]
        idata += [iidata]

    # flatten the data arrays
    odata = [e for data in odata for e in data]
    idata = [e for data in idata for e in data]

    # append all the sims per site in a single df
    dates_2_plot = []
    append = False

    for data in odata:

        if ('%s_' % (site)) in data:
            if not append:  # first occurence, declare df
                df, __ = read_csv(data)
                dates = select_dates(df, year)

                # expand to all data per day
                dates = np.repeat(dates, int(24 / (df['hod'].iloc[1] -
                                                   df['hod'].iloc[0])))[:-1]
                append = True

            else:  # append to existing df
                ndf, __ = read_csv(data)
                df = pd.concat([df, ndf], ignore_index=True)
                dates = select_dates(ndf, year)

                # expand to all data per day
                dates = np.repeat(dates, int(24 / (df['hod'].iloc[1] -
                                                   df['hod'].iloc[0])))[:-1]

            dates_2_plot += [dates]

    # masks for all years
    dates_2_plot = [e for sublist in dates_2_plot for e in sublist]

    # filter by months
    delta = df['hod'].iloc[1] - df['hod'].iloc[0]
    df.fillna(value=0., inplace=True)
    df = df.iloc[dates_2_plot]
    Dsteps = int(24. / delta)  # how many timesteps in a day?

    # exclude data 48 hours after rain (replace those gs values by 0.)
    irain = list(df[df['Ps(%s)' % ('std')].diff().gt(0) == True].index.values)

    for psi in ['psi1', 'psi2']:

        try:
            irain += list(df[df['sw(%s)' % (psi)].diff().gt(0) == True].index
                          .values)

        except KeyError:
            pass

    irain = [list(range(idx, idx + 2 * Dsteps)) for idx in np.unique(irain)]
    irain = np.unique([idx for sublist in irain for idx in sublist if
                      (idx in df.index) and (idx < len(df))])
    df.loc[irain, df.columns[df.columns.str.contains('gs')]] = 0.

    ogs = df['gs(%s)' % ('std')]
    ops = df['Ps(%s)' % ('std')]  # in MPa

    # only keep the gs 3 hours around midday
    ogs = ogs.where(np.logical_and(df['hod'] >= 9., df['hod'] <= 15.), 0.)
    ogs = running_mean(ogs, 1, delta, how='mean_non_zero')
    ops = ops.iloc[::Dsteps]  # predawn

    # make sure the lengths match
    ogs = ogs[:len(ops)]
    ops = ops[:len(ogs)]
    ops.reset_index(inplace=True, drop=True)
    ogs.reset_index(inplace=True, drop=True)

    # drop the nans where they appear!
    ops = ops[ogs.notnull()]
    ogs = ogs.dropna()

    # deal with potential outliers
    mask = ogs < ogs.quantile(.9)
    ops = ops[mask]
    ogs = ogs[mask]
    ops.reset_index(inplace=True, drop=True)
    ogs.reset_index(inplace=True, drop=True)

    # normalised gs
    ogs = ogs / np.amax(ogs.values)

    # 9th order gam fit for the gs-ps function
    gam = (LinearGAM(n_splines=11, spline_order=9, constraints='monotonic_inc')
           .gridsearch(ops.values.reshape(-1, 1), ogs.values))
    aps = np.linspace(np.amin(ops), np.amax(ops), num=500)
    ags = gam.predict(aps)
    y_int = gam.confidence_intervals(aps, width=0.95)

    ax.plot(aps, ags, linewidth=1.5, color=colour, label=site, zorder=zorder)
    axin.plot(aps, ags, linewidth=2.5, color=colour, label=site, zorder=zorder)
    ax.fill_between(aps, y_int[:, 0], y_int[:, 1], alpha=0.1, facecolor=colour,
                    zorder=zorder - 1)
    axin.fill_between(aps, y_int[:, 0], y_int[:, 1], alpha=0.1,
                      facecolor=colour, zorder=zorder - 1)
    print('min. water potential (MPa) Control at %s: %f' %
          (site, np.amin(aps)))

    for psi in ['psi1', 'psi2']:

        try:
            ogs = df['gs(%s)' % (psi)]
            ops = df['Ps(%s)' % (psi)]  # in MPa

            # only keep the gs 3 hours around midday
            ogs = ogs.where(np.logical_and(df['hod'] >= 9., df['hod'] <= 15.),
                            0.)
            ogs = running_mean(ogs, 1, delta, how='mean_non_zero')
            ops = ops.iloc[::Dsteps]  # predawn

            # make sure the lengths match
            ogs = ogs[:len(ops)]
            ops = ops[:len(ogs)]
            ops.reset_index(inplace=True, drop=True)
            ogs.reset_index(inplace=True, drop=True)

            # drop the nans where they appear!
            ops = ops[ogs.notnull()]
            ogs = ogs.dropna()

            # deal with potential outliers
            mask = ogs < ogs.quantile(.9)
            ops = ops[mask]
            ogs = ogs[mask]
            ops.reset_index(inplace=True, drop=True)
            ogs.reset_index(inplace=True, drop=True)

            # normalised gs
            ogs = ogs / np.amax(ogs.values)

            # 9th order gam fit for the gs-ps function
            gam = (LinearGAM(n_splines=11, spline_order=9,
                             constraints='monotonic_inc')
                   .gridsearch(ops.values.reshape(-1, 1), ogs.values))
            aps = np.linspace(np.amin(ops), np.amax(ops), num=500)
            ags = gam.predict(aps)
            y_int = gam.confidence_intervals(aps, width=0.95)

            ax.plot(aps, ags, linewidth=2.5, linestyle=':', color=colour,
                    zorder=-1)
            ax.fill_between(aps, y_int[:, 0], y_int[:, 1], alpha=0.1,
                            facecolor=colour, zorder=-1)
            print('min. water potential (MPa) Calib. at %s: %f' %
                  (site, np.amin(aps)))

        except KeyError:
            pass

    return


def soil_moisture_curves(figdir, sites, project1, project2=None, kmax=None,
                         plot='each', which='best'):

    """
    Retrieves the data and plots either one of the 'each' switches.

    Arguments:
    ----------
    figdir: string
        name (path) of the repository in which to output the plots

    sites: array
        names of the sites to plot

    project1: string
        name of the repository which contains the modelled data

    project2: string
        name of the repository which contains the modelled data

    kmax: array
        kmax numbers for the calibration

    plot: string
        'each' plots each individual site on a different subplot, 'all'
        stacks them all in a single plot

    which: string
        'best' picks the best performing configurations, 'worst' the
        worst ones

    Returns:
    --------
    files of the form 'what_all/each_water_stress.png' in figdir.

    """

    # colours (purple, green, orange, pink)
    col = ['#762a83', '#7fbc41', '#fdb462', '#de77ae']

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

    split_sites = split(sites, 3)
    split_projects = split(list(projects), 3)

    if kmax is not None:
        split_kmax = split(kmax, 3)

    allsites, years = find_all_combis(project1)

    if plot == 'each':  # plot each of the sites on a different axis
        plt.figure(figsize=(10., 8.))  # declare figure
        track = 0

        for i in range(3):

            Nrow = 0

            for j in range(len(split_sites[i])):

                site = split_sites[i][j]
                project = split_projects[i][j]

                if kmax is not None:
                    kn = split_kmax[i][j]

                idx = [ii for ii, e in enumerate(allsites) if e == site]
                subyears = [years[ii] for ii in idx]

                ax = plt.subplot2grid((len(split_sites[i]) * 16 + 3, 3),
                                      (Nrow, i), rowspan=16)

                if kmax is None:
                    each_water_stress(site, subyears, project, ax, col)

                else:
                    each_water_stress(site + kn, subyears, project, ax, col)

                plt.setp(ax.get_xticklabels(), visible=False)
                ax.tick_params(axis='x', length=0)

                ax.set_ylim(0., 1.)
                ax.yaxis.set_major_locator(plt.LinearLocator(3))
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

                text = ax.text(0.02, 0.9, '(%s) %s' %
                               (string.ascii_lowercase[track], site),
                               transform=ax.transAxes, weight='bold',
                               fontsize=8)
                text.set_path_effects([PathEffects.withStroke(linewidth=5,
                                                              foreground='w')])

                if i == 0:
                    ax.set_ylabel(r'g$_{\rm s, norm}$ (-)')

                if i == 0:
                    Nrow += 17

                else:
                    Nrow += 18

                track += 1

            if 'upgreek' in str(plt.rcParams['text.latex.preamble']):
                ax.set_xlabel(r'$\uptheta$ (m$^{3}$ m$^{-3}$)')

            else:
                ax.set_xlabel(r'$\theta$ (m$^{3}$ m$^{-3}$)')

        # adjust figure size
        plt.subplots_adjust(wspace=0.3, hspace=0.3)

        # add legend (only once)
        handles, labels = ax.get_legend_handles_labels()

        handles = [e for i, e in enumerate(handles) if not 'gs' in labels[i]]
        labels = [e for e in labels if not 'gs' in e]
        ax.legend(handles, labels, bbox_to_anchor=(0.225, -0.15), loc=2,
                  fancybox=True, fontsize=6, ncol=2, labelspacing=0.25)

        # save the figure
        if 'extreme' in project1:
            namefig = 'extreme_'

        if 'average' in project1:
            namefig = 'average_'

        if 'sample' in project1:
            namefig = 'calib_'

        if project2 is not None:
            namefig = '%s_' % (which)

        namefig = os.path.join(figdir, '%seach_water_stress' % (namefig))
        plt.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight')
        plt.savefig('%s.pdf' % (namefig), dpi=1800, bbox_inches='tight')
        plt.close()

    else:
        col = ['#1b9e77', '#c2a5cf', '#ffd92f', '#d95f02']

        plt.figure(figsize=(5., 4.))  # declare figure
        ax = plt.subplot(111)
        axin = inset_axes(ax, width=1.4, height=1., loc=9,
                          bbox_to_anchor=(0.105, -0.075, 1, 1),
                          bbox_transform=ax.transAxes)  # centre-top

        for line in ['top', 'bottom', 'left', 'right']:

            axin.spines[line].set_color('0.7')

        axin.tick_params(colors='0.7')

        # mimic the zoom of mark_inset function but behind the data
        ax.plot([-0.3, 0.], [0.1, 0.1], linewidth=0.9, color='0.7', zorder=-20)
        ax.plot([-0.3025, -0.54684211], [0.1015, 0.16678947], linewidth=0.9,
                color='0.7', zorder=-20)
        ax.plot([-1.231, -1.78], [0.3496, 0.5], linewidth=0.9, color='0.7',
                zorder=-20)
        ax.plot([0., -2.], [0.898, 0.81675], linewidth=0.9, color='0.7',
                zorder=-20)

        track = 0

        for i in range(len(sites)):

            site = sites[i]
            project = projects[i]

            # only plot a subset of sites
            if ((site == 'Hyytiala') or (site == 'Parco') or
               (site == 'Puechabon') or (site == 'Espirra')):

                if kmax is not None:
                    kn = kmax[i]

                idx = [ii for ii, e in enumerate(allsites) if e == site]
                subyears = [years[ii] for ii in idx]

                if kmax is None:
                    all_water_stress(site, subyears, project, ax, axin,
                                     col[track], zorder=-track)

                else:
                    all_water_stress(site + kn, subyears, project, ax, axin,
                                     col[track], zorder=-track)

                track += 1

        left, __ = ax.get_xlim()  # apply limits to the axes
        ax.set_xlim(left + 0.025, 0.)
        axin.set_xlim(-0.3, 0.)  # limit to > - 0.3 MPa
        ax.set_ylim(0., 0.9)
        axin.set_ylim(0.1, 0.895)

        ax.yaxis.set_major_locator(plt.LinearLocator(3))
        axin.yaxis.set_major_locator(plt.LinearLocator(2))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        axin.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        ax.xaxis.set_major_locator(plt.MultipleLocator(base=1.5))
        axin.xaxis.set_major_locator(plt.LinearLocator(2))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axin.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        ax.annotate('Increasingly\ndrier', xy=(0.65, 0.4),
                    xycoords='axes fraction', xytext=(0.15, 0.375),
                    ha='center', arrowprops=dict(arrowstyle="<-"))

        ax.set_ylabel(r'g$_{\rm s, norm}$ (-)', fontsize=9)
        ax.set_xlabel(r'$\Psi$$_{\rm s}$ (MPa)', fontsize=9)

        # create custom legend
        handles, labels = ax.get_legend_handles_labels()
        handles = ([Line2D([0], [0], linestyle='-', color='k', lw=0.75),
                    Line2D([0], [0], linestyle=':', color='k', lw=0.75)] +
                   handles)

        if kmax is not None:
            labels = [e.split('k')[0] for e in labels]

        labels = ['Control', r'Profit$_{\rm max}$'] + labels
        ax.legend(handles, labels, loc=2, fancybox=False, fontsize=7,
                  labelspacing=0.25)

        # save the figure
        plt.tight_layout()

        if 'extreme' in project1:
            namefig = 'extreme_'

        if 'average' in project1:
            namefig = 'average_'

        if 'sample' in project1:
            namefig = 'calib_'

        if project2 is not None:
            namefig = '%s_' % (which)

        namefig = os.path.join(figdir, '%sall_water_stress' % (namefig))
        plt.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight')
        plt.savefig('%s.pdf' % (namefig), dpi=1800, bbox_inches='tight')
        plt.close()

        return


#=======================================================================

if __name__ == "__main__":

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.weight'] = 'light'
    plt.rcParams['font.size'] = 8
    plt.rcParams['axes.titlepad'] = 0.1

    # user input
    sites = ['Hyytiala', 'Soroe', 'Loobos', 'Hesse', 'Parco', 'Puechabon',
             'Rocca1', 'Rocca2', 'ElSaler1', 'Espirra']  # sites to plot

    # define the argparse settings to read run set up file
    description = "Plot all the different site fluxes & soil \
                   variations for a given project folder"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('project1', type=str, help='project where output is')
    parser.add_argument('-c', '--project2', type=str, default=None,
                        help='project where output is')
    args = parser.parse_args()

    try:
        plt.rc('text', usetex=True)

        try:
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
                                                   r'\usepackage{upgreek}',
                                                   r'\usepackage{cmbright}']
            main(args.project1, args.project2, sites)

        except Exception as e:
            try:
                plt.rcParams['text.latex.preamble'] = \
                    [r'\usepackage{amsmath}', r'\usepackage{upgreek}']
                main(args.project1, args.project2, sites)

            except Exception as e:
                try:
                    plt.rcParams['text.latex.preamble'] = \
                        [r'\usepackage{upgreek}', r'\usepackage{cmbright}']
                    main(args.project1, args.project2, sites)

                except Exception as e:
                    try:
                        plt.rcParams['text.latex.preamble'] = \
                            [r'\usepackage{amsmath}', r'\usepackage{cmbright}']
                        main(args.project1, args.project2, sites)

                    except Exception as e:
                        try:
                            plt.rcParams['text.latex.preamble'] = \
                                [r'\usepackage{amsmath}']
                            main(args.project1, args.project2, sites)

                        except Exception as e:
                            try:
                                plt.rcParams['text.latex.preamble'] = \
                                    [r'\usepackage{upgreek}']
                                main(args.project1, args.project2, sites)

                            except Exception as e:
                                plt.rcParams['text.latex.preamble'] = \
                                    [r'\usepackage{cmbright}']
                                main(args.project1, args.project2, sites)

    except Exception as e:
        main(args.project1, args.project2, sites)
