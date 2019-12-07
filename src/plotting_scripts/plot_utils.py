# -*- coding: utf-8 -*-

"""
Functions used to support the other plotting scripts

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "useful plotting functions"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (18.01.2019)"
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
import pandas as pd  # read/write dataframes, csv files

# plotting modules
import matplotlib.dates as mdates  # format axes with dates

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM import conv  # unit converter
from TractLSM.Utils import get_main_dir  # locate files
from TractLSM.Utils import read_netcdf as netcdf  # read in data files


#=======================================================================

def find_all_combis(project, unique=True):

    """
    Finds all site-year combinations of the output.

    Arguments:
    ----------
    project: string
        name of the repository which contains the output data

    unique: boolean
        if True, combinations that are doubled (e.g. different model
        configurations for the same site-year) will be returned once
        only

    Returns:
    --------
    sites: array
        sites present in the project, in the combinations order

    years: array
        years present in the project, in the combinations order

    """

    sites = []
    years = []
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(basedir)

    # output data dir path
    oname = os.path.join(basedir, 'output')
    oname = os.path.join(os.path.join(oname, 'projects'), project)

    for file in os.listdir(oname):

        if (file.endswith('.csv')) and ('actual' in file):
            site = file.split('_')[0]
            year = int(file.split('_')[-1].split('.')[0])
            sites += [site]
            years += [year]

    if unique:
        years = [years[e] for e in range(len(sites)) if sites[e][:-1] not in
                 sites]
        sites = [e for e in sites if e[:-1] not in sites]

    return sites, years


def find_site_data(site, year, project):

    """
    Finds specific site-year files in both the input and the output for
    a project, including all model configurations.

    Arguments:
    ----------
    site: string
        individual site name to find

    year: int
        individual year to find

    project: string
        name of the repository which contains the data

    Returns:
    --------
    idata: array
        absolute paths to the input data

    odata: array
        absolute paths to the output data

    """

    idata = []
    odata = []
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    # input data dir path
    iname = os.path.join(basedir, 'input')
    iname = os.path.join(os.path.join(iname, 'projects'), project)

    # output data dir path
    oname = os.path.join(basedir, 'output')
    oname = os.path.join(os.path.join(oname, 'projects'), project)

    for file in os.listdir(iname):

        if (file.endswith('.csv')) and ('actual' in file):
            if (site in file) and (str(year) in file):
                fname = os.path.join(iname, file)
                idata += [fname]

    for file in os.listdir(oname):

        if (file.endswith('.csv')) and ('actual' in file):
            if (site in file) and (str(year) in file):
                fname = os.path.join(oname, file)
                odata += [fname]

    return idata, odata


def get_best_kmax_calib(sites, basedir):

    """
    Retrieves the best site-calibrated kmax numbers.

    Arguments:
    ----------
    sites: array
        names of the sites

    basedir: string
        name (path) of the base repository upon which the variation of
        kmax was performed

    Returns:
    --------
    kmax: array
        number names associated with the best performing config for each
        site

    """

    try:  # look for "best" kmax
        fkmax = open(os.path.join('%s_best_calib' % (basedir), 'best.txt'),
                     'r')
        flines = fkmax.readlines()
        kmax = [None] * len(sites)

        for l in flines:

            site = l.split(':')[0].strip()
            k = l.split(':')[1].strip()
            kmax[sites.index(site)] = k

        if not any(kmax):
            kmax = None

    except IOError:
        kmax = None

    return kmax


def get_best_kmax_climate(sites, basedir):

    """
    Retrieves the best type of climate for the site generation of kmax.

    Arguments:
    ----------
    sites: array
        names of the sites

    basedir: string
        name (path) of the base repository upon which the variation of
        kmax was performed

    Returns:
    --------
    climate: array
        climate type associated with the best performing site
        kmax-climate

    """

    try:  # look for "best" climate
        fclimate = open(os.path.join('%s_best_climate' % (basedir),
                        'best.txt'), 'r')
        flines = fclimate.readlines()
        climate = [None] * len(sites)

        for l in flines:

            site = l.split(':')[0].strip()
            cl = l.split(':')[1].strip()

            if 'average' in cl:
                climate[sites.index(site)] = 'Average'

            else:
                climate[sites.index(site)] = 'Extreme'

        if not any(climate):
            climate = None

    except IOError:
        climate = None

    return climate


def read_netcdf(fname, year):

    """
    Retrieves netcdf data & stores it into a dataframe

    Arguments:
    ----------
    fname: string
        input filename (with path)

    year: int
        individual year to keep

    Returns:
    --------
    df: pandas dataframe
        df containing the data

    """

    df = netcdf(fname, var_list=['Evap', 'GPP', 'LAI'])  # netcdf to df
    df = df.iloc[df.index.year == year]  # select specific year only
    df.index = range(1, len(df) + 1)

    # CABLE runs start at 12:30, so add value for 12:00
    try:
        first = [{'Evap': 0., 'GPP': 0., 'LAI': df.LAI.iloc[0]}]
        df = pd.concat([pd.DataFrame(first), df])

    except IndexError:
        try:
            first = [{'Evap': 0., 'GPP': 0., 'LAI': df.LAI.iloc[len(df)-1]}]
            df = pd.concat([pd.DataFrame(first), df])

        except IndexError:
            first = [{'Evap': 0., 'GPP': 0., 'LAI': 0.}]
            df = pd.concat([pd.DataFrame(first), df])

    return df


def select_dates(df, year):

    """
    Returns the daily indexes of the April to November subset of dates
    from any df, for a specific year.

    Arguments:
    ----------
    df: pandas dataframe
        df containing the year's surface fluxes data

    year: int
        individual year that the selection applies to

    Returns:
    --------
    The daily indexes of the subset of dates

    """

    start_date = pd.to_datetime(int(year * 1000 + df['doy'].iloc[0]),
                                format='%Y%j')

    if df['doy'].iloc[-1] > df['doy'].iloc[0]:
        end_date = pd.to_datetime(int(year * 1000 + df['doy'].iloc[-1]),
                                  format='%Y%j')

    else:
        end_date = pd.to_datetime(int((year + 1) * 1000 + df['doy'].iloc[-1]),
                                  format='%Y%j')

    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    __, iunique = np.unique(dates.month, return_index=True)
    unique_months = [dates.month[i] for i in sorted(iunique)]

    # northern hemisphere, start in April, go to November
    dates_2_plot = np.in1d(dates.month, unique_months[3:-2])

    return dates_2_plot


def convert_units(df, otype=None, var=None):

    """
    Converts data to common units across different datasets

    Arguments:
    ----------
    df: pandas dataframe
        df containing the data

    otype: string
        type of configuration looked at, e.g. std

    var: string
        in case the carbon flux is named differently to the expected
        names

    Returns:
    --------
    df: pandas dataframe
        df containing the data

    """

    if otype is None:
        try:
            df['Qle'] *= conv.Wpm2_2_mmphlfhr

        except KeyError:
            try:
                df['Evap'] *= conv.SEC_2_HLFHR

            except KeyError:
                pass

        df['GPP'] *= conv.umolCpm2ps_2_gCpm2phlfhr

        if var is not None:
            df[var] *= conv.SEC_2_HLFHR

    else:
        try:
            df['Qle'] = ((df['E(%s)' % (otype)] + df['Es(%s)' % (otype)]) *
                         conv.mmolH2Opm2ps_2_mmphlfhr)

        except KeyError:  # multicomponent case?
            df['Qle'] = (df['E_all(%s)' % (otype)] *
                         conv.mmolH2Opm2ps_2_mmphlfhr)
            df['Es(%s)' % (otype)] *= conv.mmolH2Opm2ps_2_mmphlfhr

        try:
            df['GPP'] = df['A(%s)' % (otype)] * conv.umolCpm2ps_2_gCpm2phlfhr

        except KeyError:  # multicomponent case?
            df['GPP'] = (df['A_all(%s)' % (otype)] *
                         conv.umolCpm2ps_2_gCpm2phlfhr)

        if var is not None:
            df[var] *= conv.mmolH2Opm2ps_2_mmphlfhr

    return df


def running_mean(var, N, delta, how='sum'):

    """
    Calculates a N-day running mean on either daily summed or daily
    averaged data.

    Arguments:
    ----------
    var: array
        variable on which the running mean is applied

    N: int
        length of the running mean [d]

    delta: float
        time step of the time series [h]

    how: string
        if 'sum', then cummulate the data throughout the day, otherwise
        average it

    Returns:
    --------
    var: array
        initial variable smoothed via a running mean

    """

    # How many timesteps in a day?
    Dsteps = int(24. / delta)

    if how == 'sum':
        var = var.groupby(var.index // Dsteps * Dsteps).sum()

    if how == 'mean':
        var = var.groupby(var.index // Dsteps * Dsteps).mean()

    if how == 'mean_non_zero':
        var = (var.replace(0., pd.np.nan).groupby(var.index // Dsteps * Dsteps)
                  .mean())

    if N != 1:
        var = var.rolling(N, min_periods=1).sum() / N

    return var


def arrays_boundary(var_list, which):

    """
    Calculates the upper and lower boundaries among a list of arrays,
    for each point in time

    Arguments:
    ----------
    var_list: array
        list of variables on which the min & max search is applied

    which: string
        'min' returns the lower boundary, 'max' the upper one

    Returns:
    --------
    var_boundary: array
        array comprising either all the minimum or all the maximum
        values

    """

    # boundaries for the fill between
    for i in range(1, len(var_list)):

        if i == 1:
            if which == 'min':
                var_boundary = np.minimum(var_list[0], var_list[i])

            if which == 'max':
                var_boundary = np.maximum(var_list[0], var_list[i])

        else:
            if which == 'min':
                var_boundary = np.minimum(var_boundary, var_list[i])

            if which == 'max':
                var_boundary = np.maximum(var_boundary, var_list[i])

    return var_boundary


def find_nearest(array, value):

    """
    Find the nearest value to a specific value in an array

    Arguments:
    ----------
    array: array
        list to search

    value: float
        value to search for in array

    Returns:
    --------
    Either the closest value in the array, or the array itself if the
    search fails.

    """

    array = np.asarray(array)

    try:
        idx = (np.abs(array - value)).argmin()

        return array[idx]

    except ValueError:

        return array


def split(a, N):

    """
    Splits a list or array into N-roughly equal parts.

    Arguments:
    ----------
    a: list or array
        list/array to be split, can contain any data type

    N: int
        number of sub-lists/arrays the input must be split into

    Returns:
    --------
    A list of N-roughly equal parts.

    """

    integ = int(len(a) / N)
    remain = int(len(a) % N)

    splitted = [a[i * integ + min(i, remain): (i + 1) * integ +
                min(i + 1, remain)] for i in range(N)]

    return splitted


def plot_var_2_axis(var_list, ax_list, dates=None, linestyle='-', color=None,
                    label=None, zorder=None, share_axes=False, lw=None,
                    minter=3):

    """
    Plots each of the variables in var_list on the associated axis
    ax_list.

    Arguments:
    ----------
    var_list: array
        list of variables to plot

    ax_list: array
        list of axes on which to plot

    dates: array
        list of dates associated with the data to plot

    linestyle: string
        standard linestyles

    color: string
        color used to plot the variables

    label: string
        label associated with the color and variable

    zorder: int
        zorder of the plot

    share_axes: boolean
        if True, the x axes are shared

    lw: float
        line width

    minter: int
        interval between months on the date axis [month]

    Returns:
    --------
    Plots the relevant data on the axes.

    """

    if (color == 'k'):
        if lw is None:
            lw = 1.25

    else:
        if lw is None:
            lw = 1.

    if dates is not None:

        for var in var_list:
            var.index = dates

    for i in range(len(var_list)):

        if share_axes:
            if ax_list[i] != ax_list[-1]:
                ax_list[i].get_shared_x_axes().join(ax_list[i], ax_list[-1])

        ax_list[i].plot(var_list[i].index.to_pydatetime(), var_list[i],
                        linewidth=lw, linestyle=linestyle, color=color,
                        label=label, zorder=zorder)

        if dates is not None:
            if share_axes:
                if ax_list[i] != ax_list[-1]:
                    ax_list[i].set_xticklabels([])

                else:
                    (ax_list[i].xaxis
                     .set_major_locator(mdates.MonthLocator(interval=minter)))
                    (ax_list[i].xaxis
                     .set_major_formatter(mdates.DateFormatter('%b')))

        ax_list[i].minorticks_off()

    return


def plot_var_2_var(var_list1, var_list2, ax, marker='o', ms=1.5, me=0.15,
                   color=None, label=None, zorder=None):

    """
    Plots variables as functions of other variables.

    Arguments:
    ----------
    var_list1: array
        list of variables to plot on the y axis

    var_list2: array
        list of variables to plot on the x axis

    ax: matplotlib object
        axis on which to plot

    marker: string
        marker type

    ms: float
        marker size

    me: float
        marker edge width

    color: string
        color used to plot the variables

    label: string
        label associated with the color and variable

    zorder: int
        zorder of the plot

    Returns:
    --------
    Plots the relevant data on the axis.

    """

    if zorder is None:
        zorder = -1

    ax.plot(var_list2[0], var_list1[0], marker, ms=ms, markeredgewidth=0.15,
            markeredgecolor=color, markerfacecolor='none', zorder=zorder)

    if (zorder is not None) and (zorder != -1):
        ax.plot(var_list2[1], var_list1[1], linewidth=2.5, color='k',
                zorder=zorder)
        ax.plot(var_list2[1], var_list1[1], linewidth=1.5, color=color,
                label=label, zorder=zorder)

    else:
        ax.plot(var_list2[1], var_list1[1], linewidth=2.5, color='k')
        ax.plot(var_list2[1], var_list1[1], linewidth=1.5, color=color,
                label=label)

    return


def add_bars_2_axis(var_list, ax_list, dates=None, color=None, label=None):

    """
    Adds vertical bars on the top x-axis of an axis.

    Arguments:
    ----------
    var_list: array
        list of variables to plot

    ax_list: array
        list of axes on which to plot

    dates: array
        list of dates associated with the data to plot

    color: string
        color used to plot the variables

    label: string
        label associated with the color and variable

    Returns:
    --------
    Plots the relevant data on the axes.

    """

    if dates is not None:

        for var in var_list:
            var.index = dates
            var *= -1  # so that the bars start from the top of the plot

    for i in range(len(var_list)):

            ax_list[i].bar(var.index.to_pydatetime(), var, color=color,
                           label=label, zorder=-1)

            # make the ax bars only span up to half of the vertical axis
            max_var = ax_list[i].get_ylim()[0]
            ax_list[i].set_ylim([max_var * 2., 0.])
            ax_list[i].set_yticks([np.mean(var[var < 0.]), max_var])
            real_ytick_labels = ax_list[i].get_yticks().copy()
            real_ytick_labels[real_ytick_labels < 0.] *= -1.
            ax_list[i].set_yticklabels([str(int(e)) for e in
                                       real_ytick_labels])

    return


def highlight_period(ax, period, dates):

    """
    Draws a box behind all the data plotted on an axis.

    Arguments:
    ----------
    ax: matplotlib object
        axis on which to plot

    period: array
        list of dates on which the box should span

    dates: array
        list of dates associated with the original data on the plot

    Returns:
    --------
    Plots the relevant box on the axis.

    """

    # fill background during drought, with box around
    ymin, ymax = ax.get_ylim()
    ax.fill_between(period.to_pydatetime(), ymin, ymax, color='silver',
                    facecolor='whitesmoke', linestyle='--', alpha=0.4)
    ax.set_xlim([dates.to_pydatetime()[0], dates.to_pydatetime()[-1]])
    ax.set_xticks([])

    return
