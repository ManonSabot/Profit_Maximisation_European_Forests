#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to plot the relationship of the kmax values to climate

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "plot the site-level kmax values with regard to climate"
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
from scipy import stats  # compute linear regressions

# plotting modules
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # custom legend

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.Utils import get_main_dir, read_csv  # locate & read data

# local plotting modules
from plot_utils import find_all_combis  # locate data
from plot_utils import get_best_kmax_calib  # best kmax's number name


#=======================================================================

def main(project1, project2, project3):

    """
    Main: Plots all the values of kmax across three project in
          comparison to the site-level climate.

    Arguments:
    ----------
    project1: string
        name of the repository which contains the kmax data

    project2: string
        name of the repository which contains the kmax data

    project3: string
        name of the repository which contains the kmax data

    Returns:
    --------
    'kmax_2_MAP.png' and 'kmax_LAI_2_MAP.png' in
    output/figures/final_4_paper/

    """

    # which sites are present in the data?
    sites, __ = find_all_combis(project1)
    sites = sorted(list(np.unique(np.array(sites))))

    # paths
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    basedir = os.path.join(basedir, 'output')
    figdir = os.path.join(os.path.join(basedir, 'figures'),
                          'not_shown_in_paper')

    if not os.path.isdir(figdir):  # create figdir if it doesn't exist
        os.makedirs(figdir)

    print('not comparable across sites')
    plot_kmax_climate(sites, figdir, project1, project2, project3=project3)

    figdir = os.path.join(os.path.join(basedir, 'figures'), 'final_4_paper')

    if not os.path.isdir(figdir):  # create figdir if it doesn't exist
        os.makedirs(figdir)

    print('comparable across sites')
    plot_kmax_climate(sites, figdir, project1, project2, project3=project3,
                      LAI=True)

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def kmax_values(sites, project):

    """
    Retrieves the site values of kmax present within a specific project.

    Arguments:
    ----------
    sites: array
        names of the sites to plot

    project: string
        name of the repository which contains the kmax data

    Returns:
    --------
    kmax: array
        maximum hydraulic conductance [mmol m-2 s-1 MPa-1]

    """

    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    iname = os.path.join(os.path.join(basedir, 'input'), 'projects')
    iname = os.path.join(iname, project)

    confs = ['H', '', 'L']

    if 'sample' in project:  # sites in alpha order
        confs = ['']

        # look for which is the "best" kmax
        basedir = os.path.join(os.path.join(os.path.join(basedir,
                               'output'), 'projects'), 'var_kmax')
        kmax = get_best_kmax_calib(sites, basedir)

        sites = [s + kmax[i] for i, s in enumerate(sites)]

    kmax = np.zeros((len(sites), len(confs)))

    for i in range(len(confs)):

        for j in range(len(sites)):

            for file in os.listdir(iname):

                    if (file.endswith('.csv')) and ('actual' in file) and \
                       ('%s%s_' % (sites[j], confs[i]) in file):
                        df, __ = read_csv(os.path.join(iname, file))
                        kmax[j, i] = df['kmax'].values[0]

    return kmax


def MAP_values(sites):

    """
    Retrieves the site mean annual precipication values (MAP) derived
    from the CRU datasets.

    Arguments:
    ----------
    sites: array
        names of the sites to plot

    Returns:
    --------
    MAP: array
        mean annual precipitation [mm y-1]

    """

    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    fluxsites_folder = os.path.join(os.path.join(basedir, 'input'),
                                    'fluxsites')

    info = pd.read_csv(os.path.join(fluxsites_folder, 'info_sites.csv'),
                       header=[0, 1])
    info.columns = info.columns.droplevel(level=1)
    info.index = info['Site']

    return np.array([info.loc[s, 'CRU MAP'] for s in sites])


def LAI_values(sites):

    """
    Retrieves the site composite LAI derived from the LAI climatologies
    and the CRU datasets.

    Arguments:
    ----------
    sites: array
        names of the sites to plot

    Returns:
    --------
    LAI: array
        leaf area index [m2 m-2]

    """

    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    fluxsites_folder = os.path.join(os.path.join(basedir, 'input'),
                                    'fluxsites')

    info = pd.read_csv(os.path.join(fluxsites_folder, 'info_sites.csv'),
                       header=[0, 1])
    info.columns = info.columns.droplevel(level=1)
    info.index = info['Site']

    return np.array([info.loc[s, 'Composite LAI'] for s in sites])


def plot_kmax_climate(sites, figdir, project1, project2, project3=None,
                      LAI=False):

    """
    Plots the values of kmax in comparison to either just MAP or four
    climate variables (i.e. MAP, AI, Tair, VPD).

    Arguments:
    ----------
    sites: array
        site names

    figdir: string
        name of the repository in which to output the figure

    project1: string
        name of the repository which contains the kmax data

    project2: string
        name of the repository which contains the kmax data

    project3: string
        name of the repository which contains the kmax data

    LAI: bool
        if True, plots the kmax per ground area rather than per LAI area

    Returns:
    --------
    'kmax_2_MAP.png' in figdir.

    """

    # aspect ratio
    plt.figure(figsize=(5., 4.))
    ax = plt.subplot(111)  # ax

    # plotting colors and marker shape
    msize = 200.
    alpha = 0.6

    if 'average' in project1:
        dcolours = ['#1f78b4', '#fc4e2a']

    else:
        dcolours = ['#fc4e2a', '#1f78b4']

    kmax1 = kmax_values(sites, project1)
    kmax2 = kmax_values(sites, project2)
    loops = 2

    if project3 is not None:
        dcolours = ['#7fbc41'] + dcolours
        kmax3 = kmax_values(sites, project3)
        loops = 3

    MAP = MAP_values(sites)

    if LAI:
        cLAI = LAI_values(sites)

    for i in range(loops):

        if project3 is not None:
            if i == 0:
                kmax = kmax3
                print(project3)

            elif i == 1:
                kmax = kmax1
                print(project1)

                if LAI:
                    kmax[:, 1] *= cLAI
                    kmax[:, 2] *= cLAI

            else:
                kmax = kmax2
                print(project2)

                if LAI:
                    kmax[:, 1] *= cLAI
                    kmax[:, 2] *= cLAI

        else:
            if i == 0:
                kmax = kmax1
                print(project1)

            else:
                kmax = kmax2
                print(project1)

            if LAI:
                kmax[:, 1] *= cLAI
                kmax[:, 2] *= cLAI

        if LAI:
            kmax[:, 0] *= cLAI

        if (project3 is not None) and (i == 0):
            correlate = kmax[:, 0]
            what_to = MAP

        else:
            correlate = np.concatenate((kmax[:, 0], kmax[:, 1], kmax[:, 2]))
            what_to = np.tile(MAP, 3)

        slope, intercept, r_value, p_value, __ = stats.linregress(what_to,
                                                                  correlate)

        if (i == 0) and (project3 is not None):
            size0 = [0.5 * e * msize for e in kmax[:, 0]]
            ax.scatter(MAP, kmax[:, 0], marker='o', facecolors=dcolours[i],
                       alpha=alpha, s=size0)

        else:
            idx = sites.index('Rocca2')  # avoid doubling
            size0 = [0.15 * e * msize for e in np.delete(kmax[:, 0], idx)]
            ax.scatter(np.delete(MAP, idx), np.delete(kmax[:, 0], idx),
                       marker='^', facecolors=dcolours[i], alpha=alpha,
                       s=size0)

        if (project3 is None) or ((project3 is not None) and i > 0):
            idx = sites.index('Rocca2')  # avoid doubling
            size2 = [0.15 * e * msize for e in np.delete(kmax[:, 2], idx)]
            ax.scatter(np.delete(MAP, idx), np.delete(kmax[:, 2], idx),
                       marker='v', facecolors=dcolours[i], alpha=alpha,
                       s=size2)
            size1 = [0.5 * e * msize for e in np.delete(kmax[:, 1], idx)]
            ax.scatter(np.delete(MAP, idx), np.delete(kmax[:, 1], idx),
                       marker='o', facecolors=dcolours[i], alpha=alpha,
                       s=size1)

        # highlight the sites with specific behaviours
        if ((i == 1) and (project3 is None)) or (i == 2):
            for ssite in ['Puechabon', 'Hesse']:

                idx = sites.index(ssite)

                if ssite == 'Puechabon':
                    xytext = (MAP[idx] + 25., kmax[idx, 1] + 0.02)

                if ssite == 'Hesse':
                    xytext = (MAP[idx] - 90., kmax[idx, 1] - 0.04)

                if ssite == 'Rocca1':
                    xytext = (MAP[idx] - 90., kmax[idx, 1] - 0.06)

                ax.scatter(MAP[idx], kmax[idx, 0], marker='^',
                           facecolors='none', edgecolors='k', alpha=alpha,
                           s=size0[idx])
                ax.scatter(MAP[idx], kmax[idx, 2], marker='v',
                           facecolors='none', edgecolors='k', alpha=alpha,
                           s=size2[idx])
                ax.scatter(MAP[idx], kmax[idx, 1], marker='o',
                           facecolors='none', edgecolors='k', alpha=alpha,
                           s=size1[idx])

                if ssite == 'Rocca1':
                    ax.annotate('Rocca', xy=(MAP[idx] + 100.,
                                np.amin(kmax[idx, :])), xytext=xytext)

                else:
                    ax.annotate(ssite, xy=(MAP[idx] + 100.,
                                np.amin(kmax[idx, :])), xytext=xytext)

        ax.plot([0., np.amax(MAP)],
                [intercept, intercept + np.amax(MAP) * slope],
                ':', color=dcolours[i])
        ax.set_xlim(left=300.)

        print('r: %f, r-squared: %f, and p-value: %f' % (r_value, r_value ** 2,
              p_value))

    # set N ticks on each axis
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))

    # create custom legends
    handles = [Line2D([0], [0], linestyle='', marker='^',
                      markeredgecolor='gray', fillstyle='none', markersize=6),
               Line2D([0], [0], linestyle='', marker='o',
                      markeredgecolor='gray', fillstyle='none', markersize=6),
               Line2D([0], [0], linestyle='', marker='v',
                      markeredgecolor='gray', fillstyle='none', markersize=6),
               Line2D([0], [0], linestyle=':', color=dcolours[0], marker='o',
                      alpha=alpha, markersize=6),
               Line2D([0], [0], linestyle=':', color=dcolours[1], marker='o',
                      alpha=alpha, markersize=6)]

    labels = ['High', 'Optimal', 'Low']

    if project3 is not None:
        handles += [Line2D([0], [0], linestyle=':', color=dcolours[2],
                           marker='o', alpha=alpha, markersize=6)]
        labels += ['Calibration', 'Average', 'Extreme']

    else:
        labels += ['Average', 'Extreme']

    ax.legend(handles, labels, loc=2)

    # tighten y axis
    bottom, top = ax.get_ylim()
    ax.set_ylim(0., top)

    # set legend
    ax.set_xlabel(r'{\fontsize{12pt}{3em}\selectfont{}{MAP }{\fontsize{10pt}{3em}\selectfont{}(mm y$^{-1}$)}')
    ax.set_ylabel(r'{\fontsize{12pt}{3em}\selectfont{}{k$_{\rm max}$ }{\fontsize{10pt}{3em}\selectfont{}(mmol m$^{-2}$ s$^{-1}$ MPa$^{-1}$)}')

    # save the figure
    plt.tight_layout()

    if LAI:
        namefig = os.path.join(figdir, 'kmax_per_area_2_MAP')

    else:
        namefig = os.path.join(figdir, 'kmax_per_LAI_2_MAP')

    plt.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight')
    plt.savefig('%s.pdf' % (namefig), dpi=1800, bbox_inches='tight')

    plt.close()

    return


#=======================================================================

if __name__ == "__main__":

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.weight'] = 'light'
    plt.rcParams['font.size'] = 10

    # define the argparse settings
    description = "Plot the kmax-climate relationship"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('project1', type=str,
                        help='first project where output is')
    parser.add_argument('project2', type=str,
                        help='second project where output is')
    parser.add_argument('-c', '--project3', type=str, default=None,
                        help='third project where output is')
    args = parser.parse_args()

    try:
        plt.rc('text', usetex=True)

        try:
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
                                                   r'\usepackage{cmbright}']
            main(args.project1, args.project2, args.project3)

        except Exception as e:
            try:
                plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
                main(args.project1, args.project2, args.project3)

            except Exception as e:
                plt.rcParams['text.latex.preamble'] = \
                    [r'\usepackage{cmbright}']
                main(args.project1, args.project2, args.project3)

    except Exception as e:
        main(args.project1, args.project2, args.project3)
