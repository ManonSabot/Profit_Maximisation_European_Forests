#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to plot the fluxes from the fully calibrated TractLSM

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "plot the site-level yearly calibrated fluxes at two sites"
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
import string  # automate subplot lettering

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.Utils import get_main_dir  # locate data
from TractLSM.Utils import read_csv  # read in data files
from TractLSM.SPAC import fwsoil  # plot beta functions

# local plotting modules
from plot_utils import find_site_data  # locate
from plot_utils import convert_units  # units conversions
from plot_utils import select_dates, running_mean  # transform data
from plot_utils import plot_var_2_axis  # plot


#=======================================================================

def main(project1, project2, project3, sites, years):

    """
    Main: finds the "best" performing configurations and plots them.

    Arguments:
    ----------
    project1: string
        name of repository which contains modelled data

    project2: string
        name of repository which contains modelled data

    project3: string
        name of repository which contains modelled data

    sites: array
        names of the sites to plot

    years: array
        years to plot

    Returns:
    --------
    output/figures/final_4_paper/calib_control.png

    """

    pd.options.mode.chained_assignment = None  # disable these warnings

    # paths
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    basedir = os.path.join(basedir, 'output')
    figdir = os.path.join(os.path.join(basedir, 'figures'), 'final_4_paper')
    basedir = os.path.join(basedir, 'projects')

    if not os.path.isdir(figdir):  # create figdir if it doesn't exist
        os.makedirs(figdir)

    # which calibrated variable for which project?
    projects = [project1, project2, project3]

    for project in projects:  # reorder the projects for the plots!

        if 'calib' in project:
            if 'g1' in project:
                g1 = get_best_param(sites, os.path.join(basedir, project))
                project1 = project

            elif 'fw' in project:
                fw = get_best_param(sites, os.path.join(basedir, project))
                project2 = project

        elif 'sample' in project:
            kmax = get_best_param(sites, os.path.join(basedir,
                                  'var_kmax_best_calib'))
            project3 = project

    fluxes(figdir, sites, years, project1, project2, project3, g1, fw, kmax)

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def get_best_param(sites, basedir):

    """
    Retrieves the best site-calibrated parameters.

    Arguments:
    ----------
    sites: array
        names of the sites

    basedir: string
        name (path) of the base repository upon which the variation of
        the parameter was performed

    Returns:
    --------
    params: array
        number names associated with the best performing config for each
        site

    """

    try:  # look for "best" param
        fp = open(os.path.join(basedir, 'best.txt'), 'r')
        flines = fp.readlines()
        param = [None] * len(sites)

        for l in flines:

            site = l.split(':')[0].strip()

            if site in sites:
                p = l.split(':')[1].strip()
                param[sites.index(site)] = p

        if not any(param):
            param = None

    except IOError:
        param = None

    return param


def each_flux(site, year, project, ax1, ax2, col=None, label=None):

    """
    Retrieves the data of a specific year at a single site and plots the
    fluxes.

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

    col: array string
        colour(s) in which the data must be plotted

    label: string
        data label

    Returns:
    --------
    Plots the relevant data on ax1 and ax2.

    """

    # running mean
    N = 14  # fortnightly

    # find all related data
    idata, odata = find_site_data(site, year, project)

    # plot the obs
    for data in idata:

        if ('%s_' % (site)) in data:
            df, __ = read_csv(data)  # retrieve and convert the data
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
                            share_axes=True)  # plot

    if 'list' in str(type(col)):  # if a colour array is provided
        c1 = col[0]
        c2 = col[1]

    else:
        c2 = col

    plot_ctrl = True  # makes sure the ctrl only gets plotted once

    for data in odata:

        if ('%s_' % (site)) in data:

            # retrieve and convert the data
            df, __ = read_csv(data)
            df.fillna(value=0., inplace=True)

            for j, psi in enumerate(['psi1', 'psi2']):

                try:
                    df1 = convert_units(df, otype=psi)
                    oGPP = df1['GPP']
                    oQle = df1['Qle']
                    oGPP = running_mean(oGPP, N, delta)
                    oQle = running_mean(oQle, N, delta)
                    oGPP = oGPP.iloc[dates_2_plot]
                    oQle = oQle.iloc[dates_2_plot]

                    l1 = r'Profit$_{\rm max}$(k$_{\rm calib.}$)'
                    plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                    dates=dates[dates_2_plot], color=c1,
                                    label=l1, share_axes=True)  # plot

                except KeyError:
                    pass

            if plot_ctrl:
                df1 = convert_units(df, otype='std')
                oGPP = df1['GPP']
                oQle = df1['Qle']
                oGPP = running_mean(oGPP, N, delta)
                oQle = running_mean(oQle, N, delta)
                oGPP = oGPP.iloc[dates_2_plot]
                oQle = oQle.iloc[dates_2_plot]

                l1 = r'Control %s' % (label)
                plot_var_2_axis([oGPP, oQle], [ax1, ax2],
                                dates=dates[dates_2_plot], color=c2,
                                label=l1, share_axes=True)  # plot
                plot_ctrl = False

    # tighter limits
    bottom, top = ax1.get_ylim()
    ax1.set_ylim(0., top)

    bottom, top = ax2.get_ylim()
    ax2.set_ylim(0., top)

    return


def fluxes(figdir, sites, yrs, project1, project2, project3, g1, fw, kmax):

    """
    Retrieves the data of a specific year at a single site and plots the
    fluxes.

    Arguments:
    ----------
    figdir: string
        name (path) of the repository in which to output the plots

    sites: array
        names of the sites to plot

    yrs: array
        years to plot

    project1: string
        name of repository which contains the g1 data

    project2: string
        name of repository which contains the fw data

    project3: string
        name of repository which contains the kmax data

    g1: array
        g1 numbers for the calibration

    fw: array
        fw numbers for the calibration

    kmax: array
        kmax numbers for the calibration

    Returns:
    --------
    'calib_control.png' in figdir.

    """

    # colours (green, purple, blue, orange)
    col = ['#7fbc41', '#762a83', '#4393c3', '#fdb462']

    # declare figure
    plt.figure(figsize=(4., 5.))
    rspan = 12
    cspan = 3

    # bottom axis, for beta function plots
    ax = plt.subplot2grid((rspan * 8 + 3, cspan * 5 * 3 + 10),
                          (len(sites) * 2 * (rspan + 2) + 10, cspan * 8),
                          rowspan=rspan * 2, colspan=cspan * 8)
    ax.text(0.02, 0.92, '(%s)' % (string.ascii_lowercase[2 * len(sites) *
            len(yrs)]), transform=ax.transAxes, weight='bold')

    # declare other axes (one per year at the sites)
    track = 0

    for i in range(len(sites)):

        for j in range(len(yrs)):

            if j == 0:
                ax1 = plt.subplot2grid((rspan * 8 + 3, cspan * 5 * 3 + 10),
                                       (i * 2 * (rspan + 5), j * cspan * 4),
                                       rowspan=rspan, colspan=cspan * 4)
                ax2 = plt.subplot2grid((rspan * 8 + 3, cspan * 5 * 3 + 10),
                                       (i * 2 * (rspan + 5) + rspan,
                                        j * cspan * 4), rowspan=rspan,
                                       colspan=cspan * 4, sharex=ax1)

            else:
                ax1 = plt.subplot2grid((rspan * 8 + 3, cspan * 5 * 3 + 10),
                                       (i * 2 * (rspan + 5), j * cspan * 4),
                                       rowspan=rspan, colspan=cspan * 4,
                                       sharey=ax1)
                ax2 = plt.subplot2grid((rspan * 8 + 3, cspan * 5 * 3 + 10),
                                       (i * 2 * (rspan + 5) + rspan,
                                        j * cspan * 4),
                                       rowspan=rspan, colspan=cspan * 4,
                                       sharex=ax1, sharey=ax2)

                # remove ticks that make figure look cluttered
                plt.setp(ax1.get_yticklabels(), visible=False)
                ax1.tick_params(axis='y', length=0)
                plt.setp(ax2.get_yticklabels(), visible=False)
                ax2.tick_params(axis='y', length=0)

            # remove ticks that make figure look cluttered
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.tick_params(axis='x', length=0)

            # number the subplots
            ax1.text(0.02, 0.85, '(%s)' % (string.ascii_lowercase[track]),
                     transform=ax1.transAxes, weight='bold')
            ax2.text(0.02, 0.85, '(%s)' % (string.ascii_lowercase[track+1]),
                     transform=ax2.transAxes, weight='bold')

            track += 2

            # plot the fluxes
            each_flux(sites[i] + kmax[i], int(yrs[j]), project3, ax1, ax2,
                      col=col[:2], label='ref.')
            each_flux(sites[i] + g1[i], int(yrs[j]), project1, ax1, ax2,
                      col=col[2], label=r'g$_{1 \rm calib.}$')
            each_flux(sites[i] + fw[i], int(yrs[j]), project2, ax1, ax2,
                      col=col[3], label=r'$\beta$$_{\rm calib.}$')

            if j == 0:  # number of ticks per axis
                ax1.yaxis.set_major_locator(MaxNLocator(4, integer=True,
                                                        prune='lower'))
                ax2.yaxis.set_major_locator(MaxNLocator(4, integer=True,
                                                        prune='lower'))

                # set the titles and axes names
                ax1.set_title('%s' % (sites[i]), loc='left')
                ax1.set_ylabel('GPP', fontsize=5.5)
                ax2.set_ylabel('ET', fontsize=5.5)
                ax1.get_yaxis().set_label_coords(-0.3, 0.5)
                ax2.get_yaxis().set_label_coords(-0.3, 0.5)

            # plot the beta functions
            if j == 0:
                idata, __ = find_site_data(sites[i] + fw[i], int(yrs[0]),
                                           project2)
                df, __ = read_csv(idata[0])

                p = df.iloc[0]
                sw = np.linspace(0., p.theta_sat, 500)

                if i == 0:
                    ax.plot(sw, fwsoil(p, sw), linestyle=':', color='0.5',
                            label=sites[i])

                else:
                    ax.plot(sw, fwsoil(p, sw), linestyle='--', color='0.5',
                            label=sites[i])

                ax.set_yticks([0., 0.5, 1.])
                ax.set_yticklabels(['0.0', '0.5', '1.0'])

                # set the axes names
                if 'upgreek' in str(plt.rcParams['text.latex.preamble']):
                    ax.set_xlabel(r'$\uptheta$ (m$^{3}$ m$^{-3}$)')

                else:
                    ax.set_xlabel(r'$\theta$ (m$^{3}$ m$^{-3}$)')

                ax.set_ylabel(r'$\beta$ (-)')

    # add legend
    handles, labels = ax2.get_legend_handles_labels()
    h, l = ax.get_legend_handles_labels()
    handles += h
    labels += l
    ax.legend(handles, labels, bbox_to_anchor=(-1., 0.925), loc=2,
              fancybox=True, fontsize=5.)

    # adjust figure size
    plt.subplots_adjust(wspace=0.7, hspace=1.2)

    # save the figure
    plt.savefig(os.path.join(figdir, 'calib_control.png'), dpi=1000,
                bbox_inches='tight')
    plt.savefig(os.path.join(figdir, 'calib_control.eps'), dpi=600,
                bbox_inches='tight')

    return


#=======================================================================

if __name__ == "__main__":

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.weight'] = 'light'
    plt.rcParams['font.size'] = 5.
    plt.rcParams['axes.titlepad'] = 1.8

    # user input
    sites = ['Hyytiala', 'Rocca1']  # sites to plot
    yrs = ['2002', '2003', '2005', '2006']  # years to plot

    # define the argparse settings
    description = "Plot the fluxes for several project folder"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('project1', type=str, help='project where output is')
    parser.add_argument('project2', type=str, help='project where output is')
    parser.add_argument('project3', type=str, help='project where output is')
    args = parser.parse_args()

    try:
        plt.rc('text', usetex=True)

        try:
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
                                                   r'\usepackage{cmbright}',
                                                   r'\usepackage{upgreek}']
            main(args.project1, args.project2, args.project3, sites, yrs)

        except Exception as e:
            try:
                plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
                main(args.project1, args.project2, args.project3, sites, yrs)

            except Exception as e:
                plt.rcParams['text.latex.preamble'] = \
                    [r'\usepackage{cmbright}']
                main(args.project1, args.project2, args.project3, sites, yrs)

    except Exception as e:
        main(args.project1, args.project2, args.project3, sites, yrs)
