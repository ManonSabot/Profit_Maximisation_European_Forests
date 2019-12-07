#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to plot fluxes and water stress functions for all sites
within a repository!

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "plot each year-site fluxes, gs-soil water functions"
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

try:
    from itertools import takewhile, izip  # find patterns in lists

except ImportError:
    from itertools import takewhile

import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files

# plotting modules
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator  # tick locators
import string  # automate subplot lettering

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.Utils import get_main_dir  # locate data
from TractLSM.Utils import read_csv  # read in data files

# local plotting modules
from plot_utils import find_all_combis, find_site_data  # locate
from plot_utils import convert_units  # units
from plot_utils import running_mean, arrays_boundary  # transform
from plot_utils import plot_var_2_axis, plot_var_2_var, add_bars_2_axis


#=======================================================================

def main(project):

    """
    Main: Plots each of the configurations for each of the years.

    Arguments:
    ----------
    project: string
        name of the repository which contains the modelled data

    Returns:
    --------
    files of the form 'site_fluxes_year.png' or 'site_water_stress.png'
    in output/projects/project/figures/

    """

    pd.options.mode.chained_assignment = None  # disable these warnings

    if ('sample' in project) or ('calib' in project):
        unique = False

    else:
        unique = True

    sites, years = find_all_combis(project, unique=unique)

    # paths
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    figdir = os.path.join(os.path.join(os.path.join(basedir, 'output'),
                                       'projects'), project)
    figdir = os.path.join(figdir, 'figures')

    if not os.path.isdir(figdir):  # create figdir if it doesn't exist
        os.makedirs(figdir)

    for ii in range(len(sites)):

        plot_configs_fluxes(sites[ii], years[ii], project, fill_between=True)

    for site in np.unique(sites):

        idx = [i for i, e in enumerate(sites) if e == site]
        subyears = [years[i] for i in idx]
        plot_water_stress_regs(site, subyears, project)

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def plot_configs_fluxes(site, year, project, fill_between=False):

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

    fill_between: boolean
        if True, ranges between the lowest and the highest version of
        the configuration will be plotted. This is useful for the kmax
        optimal behaviours.

    Returns:
    --------
    files of the form 'site_flux_year.png' in
    output/projects/project/figures/

    """

    # aspect ratio
    plt.figure(figsize=(5., 4.))

    # declare axes
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    ax3 = ax1.twinx()  # LAI axis
    ax4 = ax2.twinx()  # precip axis, above the drought box
    ax5 = ax1.twiny()  # drought box axis
    ax6 = ax2.twiny()  # drought box axis

    # make sure the twin axes are behind the original axes
    ax3.set_zorder(ax5.get_zorder() + 1)
    ax4.set_zorder(ax6.get_zorder() + 1)
    ax1.set_zorder(ax3.get_zorder() + 1)
    ax2.set_zorder(ax4.get_zorder() + 1)

    # number of ticks per axis
    ax1.yaxis.set_major_locator(MaxNLocator(4))
    ax2.yaxis.set_major_locator(MaxNLocator(4))
    ax3.yaxis.set_major_locator(MaxNLocator(3))
    ax5.yaxis.set_major_locator(MaxNLocator(3))

    # colours (grey, blue, purple, green, pink, orange
    col = ['dimgrey', '#4393c3', '#762a83', '#7fbc41', '#de77ae', '#fdb462']

    # running mean
    N = 14  # fortnightly

    # find all related data
    idata, odata = find_site_data(site, year, project)

    if ((any('%sH_' % (site) in e for e in idata) and
         any('%sL_' % (site) in e for e in idata)) and not
        any('%sk' % (site) in e for e in idata if map(str.isdigit, e)) and not
        any('%sRoots' % (site) in e for e in idata if map(str.isdigit, e))
            and not any('%sg1WUE' % (site) in e for e in idata if
                        map(str.isdigit, e))):
        order = ['H', '', 'L']
        labels = ['high', 'opt.', 'low']

    elif (any('%sk' % (site) in e for e in idata if map(str.isdigit, e)) or
          any('%sRoots' % (site) in e for e in idata if map(str.isdigit, e)) or
          any('%sg1WUE' % (site) in e for e in idata if map(str.isdigit, e))
          and not any('%sH_' % (site) in e for e in idata) and not
          any('%sL_' % (site) in e for e in idata)):
        cmap = plt.get_cmap('Wistia')  # update colours using oranges
        col = col[:4] + [cmap(float(e) / len(idata)) for e in
                         range(len(idata))]

        fill_between = False  # reset to false
        order = []
        labels = []

        digits = [''.join(s for s in e.split(project)[1].split('_')[0]
                  if s.isdigit()) for e in idata]
        l = [(len(d) > 1) for d in digits]

        if all(l[1:]):  # account for site names with # (e.g. Rocca1)
            digits[1:] = [e[1:] for e in digits[1:]]

        track = 0

        for digit in digits:

            try:
                if 'k%d' % (int(digit)) in idata[track]:
                    order += ['k%d' % (int(digit))]
                    labels += ['%d' % (int(digit) + 1)]

                else:
                    order += ['']
                    labels += ['ref']

            except ValueError:
                    order += ['']
                    labels += ['ref']

            track += 1

    else:
        fill_between = False  # reset to false
        order = ['']
        labels = ['opt.']

        if len(idata) == 1:
            if '%sL_' % (site) in idata[0]:
                labels = ['low']

            elif '%sH_' % (site) in idata[0]:
                labels = ['high']

            elif len(odata) > 1:
                try:
                    common_pattern = ''.join(e[0] for e in takewhile(lambda x:
                                             all(x[0] == y for y in x),
                                             izip(*odata)))

                except Exception:
                    common_pattern = ''.join(e[0] for e in takewhile(lambda x:
                                             all(x[0] == y for y in x),
                                             zip(*odata)))

                order = sorted([e.split(site)[1].split('_actual')[0] for e in
                                odata])
                labels = sorted([e.split(common_pattern)[1].split('_actual')[0]
                                 for e in odata])
                labels = [e if not 'low-' in e else e.split('low-')[1]
                          for e in labels]
                labels = [e if not e == 'low' else 'ref' for e in labels]

    if len(odata) > len(col):
        cmap = plt.get_cmap('Wistia')  # update colours using oranges
        col = col[:4] + [cmap(float(e) / len(odata)) for e in
                         range(len(odata))]

    # plot the obs
    for data in idata:

        if ('%s_' % (site)) in data:

            # retrieve and convert the data
            df, __ = read_csv(data)
            delta = df['hod'].iloc[1] - df['hod'].iloc[0]
            df[['GPP', 'Qle']] = df[['GPP', 'Qle']].where(df['PPFD'] > 0.)
            df.fillna(value=0., inplace=True)
            df = convert_units(df)

            # dates
            start_date = pd.to_datetime(year * 1000 + df['doy'].iloc[0],
                                        format='%Y%j')

            if df['doy'].iloc[-1] > df['doy'].iloc[0]:
                end_date = pd.to_datetime(year * 1000 + df['doy'].iloc[-1],
                                          format='%Y%j')

            else:
                end_date = pd.to_datetime((year + 1) * 1000 +
                                          df['doy'].iloc[-1], format='%Y%j')

            dates = pd.date_range(start=start_date, end=end_date, freq='D')
            __, iunique = np.unique(dates.month, return_index=True)
            unique_months = [dates.month[i] for i in sorted(iunique)]

            if unique_months[0] == 1:  # north. hemis., start in Feb
                dates_2_plot = np.in1d(dates.month, unique_months[1:-1])

            else:  # south. hemis., start at D14 (NaNs of running mean)
                dates_2_plot = np.in1d(dates.month, unique_months[:-1])
                dates_2_plot[:14] = False

            if unique_months[0] == 1:  # northern hemisphere
                dates_2_plot[dates.month == 1] = False  # start in Feb

            # running means on iQle and iGPP
            iGPP = running_mean(df['GPP'], N, delta)
            iQle = running_mean(df['Qle'], N, delta)
            iGPP = iGPP.iloc[dates_2_plot]
            iQle = iQle.iloc[dates_2_plot]
            plot_var_2_axis([iGPP, iQle], [ax1, ax2],
                            dates=dates[dates_2_plot], color='k', label='Obs.',
                            share_axes=True, lw=2.5, minter=2, zorder=21)

            # LAI clim with GPP on plot
            LAI = (df['LAI'].groupby(df.index // int(24. / delta) *
                                     int(24. / delta))
                            .mean()
                            .iloc[dates_2_plot])
            LAI.index = dates[dates_2_plot]
            ax3.plot(LAI.index.to_pydatetime(), LAI, color=col[0])
            ax3.tick_params(axis='y', colors=col[0])

            # daily precip from top of y axis on ET plot
            precip = (df['precip'].groupby(df.index // int(24. / delta) *
                                           int(24. / delta))
                                  .mean()
                                  .iloc[dates_2_plot])
            add_bars_2_axis([precip], [ax4], dates=dates[dates_2_plot],
                            color=col[1])
            ax4.tick_params(axis='y', colors=col[1])

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

                    plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                    dates=dates[dates_2_plot], color=col[2],
                                    label='Control', share_axes=True, lw=1.5,
                                    minter=2, zorder=22)
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
                                label = 'Profit$_{max}$'
                                plot_var_2_axis([store_GPP_psi[j][k],
                                                 store_Qle_psi[j][k]],
                                                [ax1, ax2],
                                                dates=dates[dates_2_plot],
                                                color=col[3], label=label,
                                                share_axes=True, lw=1.5,
                                                minter=2)

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
                            label = 'Profit$_{max}$(k$_{%s}$)' % (labels[i])

                            if labels == ['high', 'opt.', 'low']:
                                plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                                dates=dates[dates_2_plot],
                                                linestyle=ls, color=col[3],
                                                label=label, share_axes=True,
                                                lw=1.5, minter=2)

                            else:
                                if labels[i] == 'ref':
                                    plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                                    dates=dates[dates_2_plot],
                                                    color=col[3+i],
                                                    label=label,
                                                    share_axes=True, lw=1.5,
                                                    minter=2, zorder=23)

                                else:
                                    plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                                    dates=dates[dates_2_plot],
                                                    color=col[3+i],
                                                    label=label,
                                                    share_axes=True, lw=1.5,
                                                    minter=2)

                    except KeyError:
                        pass

    # y ticks
    ax1.locator_params(axis='y', nbins=4)
    ax2.locator_params(axis='y', nbins=4)

    # remove ticks from upper axis
    ax5.set_xticks([])
    ax6.set_xticks([])

    ax1.set_ylabel('GPP (gC.m$^{-2}$.d$^{-1}$)', fontsize=10)
    ax3.set_ylabel('LAI (m$^{2}$.m$^{-2}$)', color=col[0], fontsize=8)
    ax2.set_ylabel('ET (mm.d$^{-1}$)', fontsize=10)
    ax4.set_ylabel('Precipitation\n(mm.d$^{-1}$)', color=col[1], fontsize=8)
    ax4.yaxis.set_label_coords(1.125, 0.7)  # precip label to top y axis

    ax1.text(0.02, 0.9, '(%s)' % (string.ascii_lowercase[0]),
             transform=ax1.transAxes, weight='bold')
    ax2.text(0.02, 0.9, '(%s)' % (string.ascii_lowercase[1]),
             transform=ax2.transAxes, weight='bold')

    # add legend
    if len(idata) <= 6:
        ax2.legend(bbox_to_anchor=(1.015, -0.13), loc=3, fancybox=True,
                   fontsize=6, ncol=1, labelspacing=0.65)

    else:
        ax2.legend(bbox_to_anchor=(1.2, -0.13), loc=3, fancybox=True,
                   fontsize=6, ncol=1, labelspacing=0.65)

    plt.suptitle('%s in %d' % (site, year), fontsize=12, y=1.02)
    plt.tight_layout()

    # save the figure
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    namefig = os.path.join(os.path.join(os.path.join(basedir, 'output'),
                                        'projects'), project)
    namefig = os.path.join(namefig, 'figures')
    namefig = os.path.join(namefig, '%s_fluxes_%d' % (site, year))
    plt.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight',
                transparent=True)
    plt.savefig('%s.pdf' % (namefig), dpi=1800, bbox_inches='tight',
                transparent=True)
    plt.close()

    return


def plot_water_stress_regs(site, years, project):

    """
    Retrieves the data of for all years at a single site and plots the
    soil moisture stress there, using a polynomial fit (similar, yet not
    the same as the GAM functions in gs_vs_soil_moisture).
    Arguments:
    ----------
    site: string
        individual site name
    years: array
        all data years to plot
    project: string
        name of the repository which contains the data to plot
    Returns:
    --------
    files of the form 'site_water_stress.png' in
    output/projects/project/figures/
    """

    # aspect ratio
    plt.figure(figsize=(5., 4.))

    # axes
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    # colours (purple, green, pink, orange
    col = ['#762a83', '#7fbc41', '#de77ae', '#fdb462']
    mk = 'o'

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

    if ((any('%sH' % (site) in e for e in idata) and
       any('%sL' % (site) in e for e in idata)) and not
       any('%sk%d' % (site, str.isdigit(e)) in e for e in idata) and not
       any('%sRoots%d' % (site, str.isdigit(e)) in e for e in idata) and not
       any('%sg1WUE%d' % (site, str.isdigit(e)) in e for e in idata)):
        order = ['H', '', 'L']
        labels = ['high', 'opt.', 'low']

        # reverse green and pink to keep green consistency for kmax
        green = col[1]
        pink = col[2]
        col[1] = pink
        col[2] = green

    elif (any('%sk%d' % (site, str.isdigit(e)) in e for e in idata) or
          any('%sRoots%d' % (site, str.isdigit(e)) in e for e in idata) or
          any('%sg1WUE%d' % (site, str.isdigit(e)) in e for e in idata) and not
          (any('%sH_' % (site) in e for e in idata) or
           any('%sL_' % (site) in e for e in idata))):
        cmap = plt.get_cmap('Wistia')  # update colours using oranges
        col = col[:2] + [cmap(float(e) / len(idata)) for e in
                         range(len(idata))]
        order = []
        labels = []

        for e in idata:

            digit = ''.join(s for s in e.split(project)[1].split('_')[0]
                            if s.isdigit())

            if len(digit) > 1:
                digit = digit[1]  # account for site names with number

            try:
                if 'k%d' % (int(digit)) in e:
                    order += ['k%d' % (int(digit))]
                    labels += ['%d' % (int(digit) + 1)]

                else:
                    order += ['']
                    labels += ['ref']

            except ValueError:
                    order += ['']
                    labels += ['ref']

    elif (any('%sk%d' % (site, str.isdigit(e)) in e for e in idata) and
          (any('%sH' % (site) in e for e in idata) or
           any('%sL' % (site) in e for e in idata))):
        order = ['', ] * len(idata)
        labels = ['', ] * len(idata)

    else:
        order = ['']
        labels = ['opt.']

    if len(odata) > len(col):
        cmap = plt.get_cmap('Wistia')  # update colours using oranges
        col = col[:4] + [cmap(float(e) / len(odata)) for e in
                         range(len(odata))]

    # append all the sims per config in a single df
    dfs = []

    for i in range(len(order)):

        append = False

        for data in odata:

            conf = order[i]

            if ('%s%s_' % (site, conf)) in data:
                if not append:  # first occurence, declare df
                    df, __ = read_csv(data)

                else:  # append to existing df
                    ndf, __ = read_csv(data)
                    df = pd.concat([df, ndf], ignore_index=True)

                append = True

        dfs += [df]

    # plot gs from the sim(s)
    for i in range(len(order)):

        df = dfs[i]
        conf = order[i]

        # retrieve and convert the data
        delta = df['hod'].iloc[1] - df['hod'].iloc[0]
        df.fillna(value=0., inplace=True)

        Dsteps = int(24. / delta)  # how many timesteps in a day?

        if i == 0:
            ogs = df['gs(%s)' % ('std')]
            osw = df['sw(%s)' % ('std')]

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
            mask = osw < osw.quantile(.95)
            osw = osw[mask]
            ogs = ogs[mask]
            osw.reset_index(inplace=True, drop=True)
            ogs.reset_index(inplace=True, drop=True)

            # normalised gs
            ngs = ogs / np.amax(ogs)

            # 5th order polynomial fit for the gs-sw function
            f = np.poly1d(np.polyfit(osw, ogs, 3))
            asw = np.linspace(np.amin(osw), np.amax(osw), len(osw) * 4)
            ags = f(asw)

            # 5th order polynomial fit for the ngs-sw function
            f = np.poly1d(np.polyfit(osw, ngs, 3))
            asw = np.linspace(np.amin(osw), np.amax(osw), len(osw) * 4)
            angs = f(asw)

            plot_var_2_var([ogs, ags], [osw, asw], ax1, ms=3., me=0.5,
                           color=col[0], zorder=21)  # plot
            plot_var_2_var([ngs, angs], [osw, asw], ax2, ms=3., me=0.5,
                           color=col[0], label='Control', zorder=21)

        for psi in ['psi1', 'psi2']:

            try:
                ogs = df['gs(%s)' % (psi)]
                osw = df['sw(%s)' % (psi)]

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
                mask = osw < osw.quantile(.95)
                osw = osw[mask]
                ogs = ogs[mask]
                osw.reset_index(inplace=True, drop=True)
                ogs.reset_index(inplace=True, drop=True)

                # normalised gs
                ngs = ogs / np.amax(ogs)

                # 3rd order polynomial fit for the gs-sw function
                f = np.poly1d(np.polyfit(osw, ogs, 3))
                asw = np.linspace(np.amin(osw), np.amax(osw), len(osw) * 4)
                ags = f(asw)

                # 3rd order polynomial fit for the ngs-sw function
                f = np.poly1d(np.polyfit(osw, ngs, 3))
                asw = np.linspace(np.amin(osw), np.amax(osw), len(osw) * 4)
                angs = f(asw)

                # plot
                label = 'Profit$_{max}$(k$_{%s}$)' % (labels[i])

                if labels[i] == 'ref':
                    plot_var_2_var([ogs, ags], [osw, asw], ax1, marker=mk,
                                   ms=3., me=0.5, color=col[1+i], zorder=22)
                    plot_var_2_var([ngs, angs], [osw, asw], ax2, marker=mk,
                                   ms=3., me=0.5, color=col[1+i], label=label,
                                   zorder=22)

                else:
                    plot_var_2_var([ogs, ags], [osw, asw], ax1, marker=mk,
                                   ms=3., me=0.5, color=col[1+i])
                    plot_var_2_var([ngs, angs], [osw, asw], ax2, marker=mk,
                                   ms=3., me=0.5, color=col[1+i], label=label)

            except KeyError:
                pass

    ax1.set_ylabel('gs (mol.m$^{-2}$.s$^{-1}$)', fontsize=10)
    ax2.set_ylabel('gs$_{norm}$ (-)', fontsize=10)
    ax2.set_xlabel(r'$\theta$ (m$^{3}$.m$^{-3}$)', fontsize=10)

    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_xticklabels([])

    ax1.text(0.02, 0.9, '(%s)' % (string.ascii_lowercase[0]),
             transform=ax1.transAxes, weight='bold')
    ax2.text(0.02, 0.9, '(%s)' % (string.ascii_lowercase[1]),
             transform=ax2.transAxes, weight='bold')

    # add legend
    if len(idata) <= 6:
        ax2.legend(bbox_to_anchor=(1.015, 0.), loc=3, fancybox=True,
                   fontsize=6, ncol=1, labelspacing=0.65)

    plt.suptitle('%s' % (site), fontsize=12, y=1.02)
    plt.tight_layout()
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    namefig = os.path.join(os.path.join(os.path.join(basedir, 'output'),
                                        'projects'), project)
    namefig = os.path.join(namefig, 'figures')
    namefig = os.path.join(namefig, '%s_water_stress' % (site))
    plt.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight',
                transparent=True)
    plt.savefig('%s.pdf' % (namefig), dpi=1800, bbox_inches='tight',
                transparent=True)
    plt.close()

    return


#=======================================================================

if __name__ == "__main__":

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.weight'] = 'light'
    plt.rcParams['font.size'] = 9
    plt.rc('text', usetex=True)

    # define the argparse settings to read run set up file
    description = "Plot all the different site fluxes & soil \
                   variations for a given project folder"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('project', type=str,
                        help='project where the output is')
    args = parser.parse_args()

main(args.project)
