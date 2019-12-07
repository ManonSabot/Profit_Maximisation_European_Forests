#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to plot the ProfitMax algorith behaviour on both
instantaneous and longer timescales. This is a useful tool which can be
played with to better understand model behaviour.

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "explore the ProfitMax algorithm behaviour"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (08.10.2018)"
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
import matplotlib.pyplot as plt
import string  # automate subplot lettering

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM import dparams  # model parameters
from TractLSM.SPAC import hydraulics, Weibull_params  # hydraulics
from TractLSM.SPAC import net_radiation  # big leaf canopy radiation
from TractLSM.CH2OCoupler.ProfitMax import hydraulic_cost, photo_gain
from TractLSM.TraitCoordination import optimal_kmax  # long term


#=======================================================================

def main():

    """
    Main: Explores the ProfitMax algorith behaviour on both
          instantaneous and longer timescales of optimisation. Plots and
          saves outputs.

    Returns:
    --------
    files in output/figures/final_4_paper/ and in
    output/figures/not_shown_in_paper/

    """

    # how the model works
    kmax_model_behaviour()

    # model's sensitivity to Pcrit
    rcrit_model_behaviour(how='both')

    # kmax dependency on climate
    kmax_2_var(var='both')

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def declared_params():

    """
    Sets up the model to run, by making the parameter class a pandas
    series and adding missing forcings to it.

    Returns:
    --------
    p: pandas series
        met forcings & params necessary to run the model

    """

    # make the default param class a pandas series
    p = dparams
    attrs = vars(p)
    p = {item[0]: item[1] for item in attrs.items()}
    p = pd.Series(p)

    # add missing params, radiation is arbritrary (model not sensitive)
    p = p.append(pd.Series([25., 1., 1., p.albedo_ws],
                           index=['Tair', 'VPD', 'scale2can', 'albedo_s']))
    p = p.append(pd.Series([1000.], index=['PPFD']))
    p = p.append(pd.Series([net_radiation(p)], index=['Rnet']))

    return p


def opt_stream(p, res, photo, case, parallel, solstep, symbolic):

    """
    Finds the instateneous profit maximization, following the
    optmization criterion for which, at each instant in time, the
    stomata regulate canopy gas exchange and pressure to achieve the
    maximum profit, which is the maximum difference between the
    normalized photosynthetic gain (gain) and the hydraulic cost
    function (cost). That is when d(gain)/dP = d(cost)/dP.

    Arguments:
    ----------
    p: pandas series
        met forcings & params

    res: string
        either 'low' (default), 'med', or 'high' to run the optimising
        solver

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    case: int
        either 1, where Ci isn't stricktly bound by the standard model,
        or 2 where it is

    parallel: boolean or int
        True parallelises the minimization solving of Ci over the
        system's CPUs, an int input will parallelise over N=int CPUs

    solstep: string
        'var' means the stepping increment between the lower Ci and Cs
        changes depending on Cs (e.g. there are N fixed values between
        the lower Ci and Cs), while 'fixed' means the stepping increment
        between the lower Ci and Cs is fixed (e.g. N values between the
        lower Ci and Cs is not fixed)s

    symbolic: boolean
        if True, a symbolic expression of Ci will be used to solve for
        it

    Returns:
    --------
    P: array
        leaf water potential [MPa], from the soil water potential Ps to
        the critical water potential Pcrit for which cavitation of the
        xylem occurs

    cost: array
        hydraulic cost [unitless], depending on the plant vulnerability
        curve

    gain: array
        unitless instantaneous photosynthetic gains for possible values
        of Ci minimized over the hydraulic stream

    """

    # hydraulics
    P, trans = hydraulics(p, res=res)
    cost, __ = hydraulic_cost(p, P)

    # assimilation
    gain, __ = photo_gain(p, trans, photo, res, case, parallel, solstep,
                          symbolic)

    return P, cost, gain


def Px(p, x):

    """
    Finds the leaf water potential associated with a specific x%
    decrease in hydraulic conductance, using the plant vulnerability
    curve.

    Arguments:
    ----------
    p: pandas series
        met forcings & params

    x: float
        percentage loss in hydraulic conductance

    Returns:
    --------
    Px: float
        leaf water potential [MPa] at which x% decrease in hydraulic
        conductance is observed

    """

    b, c = Weibull_params(p)  # MPa, unitless
    Px = -b * ((- np.log(1 - float(x) / 100.)) ** (1. / c))

    return Px


def plot_model(ax, P, gain, cost, colours, label=None, alpha=None,
               linestyle='-', zorder=None):

    """
    Plots each of the gain, cost, and net profit functions, along the
    water potential stream.

    Arguments:
    ----------
    ax: matplotlib object
        axis on which to plot

    P: array
        leaf water potential [MPa], from the soil water potential Ps to
        the critical water potential Pcrit for which cavitation of the
        xylem occurs

    gain: array
        unitless instantaneous photosynthetic gains for possible values
        of Ci minimized over the hydraulic stream

    cost: array
        hydraulic cost [unitless], depending on the plant vulnerability
        curve

    colours: array
        color used to plot gain, cost, and profit

    label: string
        label associated with the specific run mode

    alpha: float
        transparency

    linestyle: string
        standard linestyles

    zorder: int
        zorder of the plot

    Returns:
    --------
    Plots the relevant data on the axis.

    """

    ax.plot(-P, gain, linestyle, color=colours[2], label=label, alpha=alpha,
            linewidth=1.5, zorder=zorder)
    ax.plot(-P, cost, linestyle, color=colours[0], alpha=alpha,
            linewidth=1.5, zorder=zorder)
    ax.plot(-P, gain - cost, linestyle, color=colours[1], alpha=alpha,
            linewidth=2.5, zorder=zorder)

    return


def behavioural_markers(ax, P, profits, costs):

    """
    Plots the maximum profit marks.

    Arguments:
    ----------
    ax: matplotlib object
        axis on which to plot

    P: array
        leaf water potential [MPa], from the soil water potential Ps to
        the critical water potential Pcrit for which cavitation of the
        xylem occurs

    profits: array
        net profit, i.e. gain - cost for all three strategies

    costs: array
        hydraulic cost [unitless], depending on the plant vulnerability
        curve, for all three strategies

    Returns:
    --------
    Highlights the profit marks.

    """

    # most profit?
    Popts = [-P[np.argmax(profits[0])], -P[np.argmax(profits[-1])]]

    # plot the arrow
    ax.annotate('', xy=(Popts[-1], np.amax(profits[-1]) - 0.015),
                xytext=(Popts[0], np.amax(profits[0]) + 0.015),
                xycoords='data', textcoords='data',
                arrowprops=dict(arrowstyle='<->'), fontsize=10., ha='right',
                va='bottom')

    # write the text
    ax.text(Popts[-1], np.amax(profits[0]) + 0.02,
            r'$\Delta$$\Psi$$_{\rm leaf,opt}$ = ' +
            r'{0:.1f} MPa'.format(Popts[-1] - Popts[0]), fontsize=9.,
            ha='left', va='center')

    return


def annotate_model(ax, colours):

    """
    Annotates the modelled curves.

    Arguments:
    ----------
    ax: matplotlib object
        axis on which to plot

    colours: array
        color used to plot gain, cost, and profit

    Returns:
    --------
    Plots the relevant data on the axis.

    """

    ax.text(0.7, 0.9, 'Carbon\ngain', color=colours[2], fontsize=12.,
            va='center', ha='center', transform=ax.transAxes)
    ax.text(0.9, 0.8, 'Hydraulic\ncost', color=colours[0], fontsize=12.,
            va='center', ha='center', transform=ax.transAxes)
    ax.text(0.75, 0.315, 'Net\nprofit', color=colours[1], fontsize=12.,
            va='center', ha='center', transform=ax.transAxes)

    return


def Tair_from_VPD(vpd, rh):

    """
    Calculates the air temperature that can physically be associated
    with specific values of vapour pressure deficit, within a range of
    possible relative humidity.

    Arguments:
    ----------
    vpd: array or float
        vapour pressure deficit [kPa]

    rh: array or float
        relative humidity [unitless]

    Returns:
    --------
    Tair: array
        near surface air temperature [deg C]

    """

    rh *= 0.01  # from % to 0.-1.9
    Tair = (-237.3 * np.log(vpd / (0.61078 * (1. - rh))) /
            (np.log(vpd / (0.61078 * (1. - rh))) - 17.27))  # Tetens

    return Tair


def stream_ticks(p, P, Pcrit=False):

    """
    Sets up the ticks and tick labels on the x axis, following commonly
    referred to water potentials on the hydraulic stream.

    Arguments:
    ----------
    p: pandas series
        met forcings & params

    P: array
        leaf water potential [MPa], from the soil water potential Ps to
        the critical water potential Pcrit for which cavitation of the
        xylem occurs

    Pcrit: boolean
        if True, the 'Pcrit' is added to the tick list

    Returns:
    --------
    Renders the x axis accordingly.

    """

    P12 = Px(p, 12)
    iP12 = np.argmin(P - P12 >= 0.)

    Pticks = [-np.amax(P), -P[iP12], p.P50, p.P88]
    Ptick_labels = [r'$\Psi$$_{\rm sat}$', r'$\Psi$$_{\rm 12}$',
                    r'$\Psi$$_{\rm 50}$', r'$\Psi$$_{\rm 88}$']

    if Pcrit:
        Pticks += [-P[-1]]
        Ptick_labels += [r'$\Psi$$_{\rm crit}$']

    return Pticks, Ptick_labels


def format_model_axes(ax):

    """
    Formats the x and y axes, for better rendering.

    Arguments:
    ----------
    ax: matplotlib object
        axis on which to plot

    Returns:
    --------
    Renders the plot accordingly.

    """

    # make sure both x and y axes start at 0.
    ax.autoscale(enable=True, axis='both', tight=True)
    ax.set_xlabel(r'$\Psi$ (-MPa)', fontsize=12.)
    ax.set_ylabel(r'Gain ${\vert}$ Cost ${\vert}$ Profit', fontsize=14.)

    return


def get_fig_dir():

    """
    Returns the figure directory in which to store the plots

    """

    basedir = os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))

    while 'src' in basedir:
        basedir = os.path.dirname(basedir)

    fig_dir = os.path.join(os.path.join(basedir, 'output'), 'figures')

    if not os.path.isdir(os.path.dirname(fig_dir)):
        os.makedirs(os.path.dirname(fig_dir))

    return fig_dir


def kmax_model_behaviour():

    """
    Explores the ProfitMax algorith behaviour on the instantaneous
    timescale of optimisation, following three potential ways in which
    hydraulic and photosynthetic traits could be coordinated on longer
    time scales.

    Returns:
    --------
    output/figures/final_4_paper/profitmax_behaviours.png

    """

    # figure size and ratio
    fig = plt.figure(figsize=(5., 4.))
    ax = plt.subplot(111)

    # retrieve the three potential "behavioural" kmax values
    p = declared_params()
    kmax1, kmax2, kmax3 = optimal_kmax(p, photo, strategy='all')
    p = declared_params()

    p.kmax = kmax1[0]  # high behaviour?
    P, cost_high, gain_high = opt_stream(p, res, photo, case, parallel,
                                         solstep, symbolic)
    p = declared_params()

    p.kmax = kmax3[0]  # low behaviour?
    P, cost_low, gain_low = opt_stream(p, res, photo, case, parallel,
                                       solstep, symbolic)
    p = declared_params()

    p.kmax = kmax2[0]  # optimal behaviour?
    P, cost, gain = opt_stream(p, res, photo, case, parallel, solstep,
                               symbolic)

    lines = ['-', '--', ':']  # line types for the plots
    labels = [r'k$_{\rm max,opt}$', r'k$_{\rm max,high}$',
              r'k$_{\rm max,low}$']

    # plot the three alternative instantaneous optimisation
    plot_model(ax, P, gain_high, cost_high, colours, linestyle=lines[1],
               label=labels[1])
    plot_model(ax, P, gain, cost, colours, linestyle=lines[0], label=labels[0])
    plot_model(ax, P, gain_low, cost_low, colours, linestyle=lines[2],
               label=labels[2])

    # annotate the change in Pleaf,opt with an arrow
    behavioural_markers(ax, P, [gain_high - cost_high, gain - cost,
                        gain_low - cost_low], [cost_high, cost, cost_low])

    # annotate the model curves
    annotate_model(ax, colours)

    # retrieve horizontal axes ticks and labels
    Pticks, Ptick_labels = stream_ticks(p, P, Pcrit=True)

    # set axes ticks and labels
    ax.set_xticks(Pticks)
    ax.set_xticklabels(Ptick_labels)
    ax.set_yticks([0., 0.5, 1.])  # cost|gain|profit ticks
    ax.set_yticklabels(['0.0', '0.5', '1.0'])

    # make sure both x and y axes start at 0.
    format_model_axes(ax)
    plt.tight_layout()

    # set the legend and modify it
    lgd = ax.legend(fontsize=10., frameon=False, numpoints=2., borderpad=0.,
                    handletextpad=0.2, loc=2, bbox_to_anchor=(-0.01, 0.995))

    for handle in lgd.legendHandles:
        handle.set_linewidth(1.)
        handle.set_color('dimgrey')

    for label in lgd.get_texts():
        label.set_color('dimgrey')

    # save the figure
    figdir = os.path.join(get_fig_dir(), 'final_4_paper')
    namefig = os.path.join(figdir, 'profitmax_behaviours')
    fig.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight')
    fig.savefig('%s.eps' % (namefig), dpi=600, bbox_inches='tight')
    plt.close()

    return


def rcrit_model_behaviour(how=None):

    """
    Explores the ProfitMax algorith behaviour on the instantaneous
    timescale of optimisation, by changing the stomate's control on the
    optimisation (through Pcrit, via ratio_crit). Hydraulic and
    photosynthetic traits are assumed to be coordinated on longer time
    scales.

    Arguments:
    ----------
    how: string
        if None the three variations on stomatal control are plotted
        without any adjustment of kmax depending on ratio_crit, if
        'kmax' then kmax is adjusted for each of these ratio_crit
        values. If 'both', all of the above.

    Returns:
    --------
    files are: output/figures/ratiocrits.png,
               output/figures/kmax_ratiocrits.png,
               output/figures/final_4_paper/profitmax_ratiocrits.png.

    """

    # figure size and ratio
    fig = plt.figure(figsize=(5., 4.))
    ax = plt.subplot(111)
    axes = [ax]

    if how == 'both':
        fig = plt.figure(figsize=(10., 4.25))
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122, sharey=ax1)
        plt.subplots_adjust(wspace=0.125)
        axes = [ax1, ax2]

    # declare three variations on stomatal control
    p = declared_params()
    rcrits = [p.ratiocrit, 0.15, 0.45]
    lines = ['-', '--', ':']  # line types for the plots

    # the variations on stomatal control
    track = 0

    for ax in axes:  # loop over the two sub-plots

        Popt = []  # empty array to store Popt values

        for i in range(len(rcrits)):  # loop over the three variations

            p = declared_params()  # reset the parameters each time
            p.ratiocrit = rcrits[i]
            label = 'r$_{crit}$' + '={0:.2f}'.format(rcrits[i])

            if i == 0:  # get optimal kmax
                __, kmax, __ = optimal_kmax(p, photo, strategy='optimal')
                p.kmax = kmax[0]

            if (how == 'kopt') or (how == 'both' and track == 1):
                if i > 0:  # if not default control
                    __, kmax, __ = optimal_kmax(p, photo, strategy='optimal')
                    p.kmax = kmax[0]

            # behaviour?
            P, cost, gain = opt_stream(p, res, photo, case, parallel, solstep,
                                       symbolic)
            plot_model(ax, P, gain, cost, colours, label=label,
                       linestyle=lines[i])

            # most profit?
            iopt = np.argmax(gain - cost)
            ax.text(-P[iopt], np.amax(gain - cost) * 1.0025, str(i + 1),
                    color='k', fontsize=10., va='center', ha='center',
                    zorder=10)
            Popt += [-P[iopt]]

            if (how == 'both') and (i == 0):  # letter the subplots
                ax.text(0.02, 0.95, '(%s)' % (string.ascii_lowercase[track]),
                        transform=ax.transAxes, fontweight='heavy',
                        fontsize=12)

            if i == 2:  # set the ticks and labels
                Pticks, Ptick_labels = stream_ticks(p, P)
                ax.set_xticks(Pticks)
                ax.set_xticklabels(Ptick_labels)
                ax.set_yticks([0., 0.5, 1.])
                ax.set_yticklabels(['0.0', '0.5', '1.0'])

        # write the Pleaf,opt on the plot
        info = r'1. P$_{\rm leaf,opt}$' + ' = -{0:.2f} MPa'.format(Popt[0]) + \
               '\n' + \
               r'2. P$_{\rm leaf,opt}$' + ' = -{0:.2f} MPa'.format(Popt[1]) + \
               '\n' + \
               r'3. P$_{\rm leaf,opt}$' + ' = -{0:.2f} MPa'.format(Popt[2])
        ax.text(0.7, 0.415, info, color='k', fontsize=8., va='bottom',
                ha='left', transform=ax.transAxes)  # long string

        # make sure both axes start at 0.
        format_model_axes(ax)

        # next iteration on loop
        track += 1

    if how == 'both':
        ax2.set_ylabel('')  # do not repeat cost|gain|profit ticks
        ax = ax1

    # set the legend and modify it
    lgd = ax.legend(fontsize=10., frameon=False, numpoints=2., borderpad=0.,
                    handletextpad=0.2, loc=2, bbox_to_anchor=(-0.15, -0.075))

    for handle in lgd.legendHandles:
        handle.set_linewidth(1.)
        handle.set_color('dimgrey')

    for label in lgd.get_texts():
        label.set_color('dimgrey')

    # save the figure
    if how is None:
        namefig = os.path.join(get_fig_dir(), 'ratiocrits')

    if how is not None:
        if how == 'kopt':
            namefig = os.path.join(get_fig_dir(), 'kmax_ratiocrits')

        else:
            plt.tight_layout()
            figdir = os.path.join(get_fig_dir(), 'not_shown_in_paper')

            namefig = os.path.join(figdir, 'profitmax_ratiocrits')

    fig.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight')
    fig.savefig('%s.eps' % (namefig), dpi=600, bbox_inches='tight')
    plt.close()

    return


def kmax_2_var(var='Vcmax'):

    """
    Explores the coordination between hydraulic and photosynthetic
    traits assumed on longer timescales of optimisation, by varying
    drivers of the assumed coordination.

    Arguments:
    ----------
    var: string
        defines the first variable to vary (e.g. Vcmax, Tair, CO2,
        etc.). If 'both' then the relationship to both Vcmax and Tair is
        returned.

    Returns:
    --------
    files are: output/figures/kmax_2_var.png,
               output/figures/final_4_paper/profitmax_2_forcings.png.

    """

    # figure size and ratio
    fig = plt.figure(figsize=(5., 4.))
    ax = plt.subplot(111)
    axes = [ax]
    variables = [var]

    if var == 'both':
        fig = plt.figure(figsize=(8., 3.5))
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)
        plt.subplots_adjust(wspace=0.125)
        axes = [ax1, ax2]
        variables = ['Vcmax', 'Tair']

    # declare VPD options to loop on
    VPD = np.array([0.5, 1., 2., 3.])
    ls = ['-', '--', '-.', ':']

    for i in range(len(variables)):  # loop on the tested variables

        p = declared_params()  # reset the parameters
        ivar = variables[i]
        ax = axes[i]

        if ivar == 'Vcmax':
            Vcmax25 = list(np.arange(20., 151., 10.))
            Trange = [25.]

            kmax = np.ones((len(VPD), len(Vcmax25)))

        elif ivar == 'Tair':
            Trange = np.arange(5., 41.)
            Vcmax25 = [p.Vmax25]

            # bound the temperature for it to be valid
            RH = [5., 95.]
            Tair = np.ones((len(RH), len(VPD)))

            for i in range(len(RH)):
                Tair[i, :] = Tair_from_VPD(VPD, RH[i])

            kmax = np.ones((len(VPD), len(Trange)))

        for j in range(len(Trange)):  # restrict valid temperatures

            p = declared_params()  # reset the parameters
            p.Tair = Trange[j]
            __, kmax_opt, __ = optimal_kmax(p, photo, strategy='optimal',
                                            Vcmax25=Vcmax25, VPD=list(VPD))

            for k in range(len(VPD)):

                if ivar == 'Vcmax':

                    for l in range(len(Vcmax25)):

                        kmax[k, l] = kmax_opt[l, k, 0]

                elif ivar == 'Tair':
                    kmax[k, j] = kmax_opt[k, 0]

        for j in range(len(VPD)):  # loop on the VPD values

            kmax[j, :][kmax[j, :] <= 1.e-9] = np.nan

            if len(np.argwhere(np.isnan(kmax[j, :]))) > 0:
                jmax = np.argwhere(np.isnan(kmax[j, :]))[-1][0]
                kmax[j, jmax:] = np.nan

            if ivar == 'Vcmax':
                ax.plot(Vcmax25, kmax[j, :] * p.LAI, color='k',
                        linestyle=ls[j], linewidth=1.5,
                        label='D = %s kPa' % str(VPD[j]))
                ax.set_xlabel(r'{\fontsize{12pt}{3em}\selectfont{}{V$_{\rm cmax,25}$ }{\fontsize{10pt}{3em}\selectfont{}($\mu$mol m$^{-2}$ s$^{-1}$ )}')

            elif ivar == 'Tair':
                if Tair[0, j] > Trange[0]:
                    jmin = np.argmax(Tair[0, j] <= Trange)
                    kmax[j, :jmin] = np.nan

                if Tair[1, j] < Trange[-1]:
                    jmax = np.argmin(Tair[1, j] >= Trange)
                    kmax[j, jmax:] = np.nan

                ax.plot(Trange, kmax[j, :] * p.LAI, color='k', linestyle=ls[j],
                        linewidth=1.5, label='D = %s kPa' % str(VPD[j]))
                ax.set_xlabel(r'{\fontsize{12pt}{3em}\selectfont{}{T$_{\rm air}$ }{\fontsize{10pt}{3em}\selectfont{}($^\circ$C)}')

            if var == 'both':  # letter the subplots
                ax.text(0.02, 0.95, '(%s)' % (string.ascii_lowercase[i]),
                        transform=ax.transAxes, fontweight='heavy',
                        fontsize=10)

        # reduce the number of ticks on each axis
        ax.locator_params(axis='x', nbins=5)
        ax.locator_params(axis='y', nbins=4)

        # make sure both axes start at 0.
        ax.autoscale(enable=True, axis='both', tight=True)

    # set the legend and modify it
    if var == 'both':
        ax1.set_ylabel(r'{\fontsize{12pt}{3em}\selectfont{}{k$_{\rm max,opt}$ }{\fontsize{10pt}{3em}\selectfont{}(mmol m$^{-2}$ s$^{-1}$ MPa$^{-1}$)}')
        lgd = ax1.legend(fontsize=8., frameon=False, numpoints=0.8,
                         borderpad=0., handletextpad=0.4, loc=2,
                         bbox_to_anchor=(0.01, 0.9375))

    else:
        ax.set_ylabel(r'{\fontsize{12pt}{3em}\selectfont{}{k$_{\rm max,opt}$ }{\fontsize{10pt}{3em}\selectfont{}(mmol m$^{-2}$ s$^{-1}$ MPa$^{-1}$)}')
        lgd = ax.legend(fontsize=7., frameon=False, numpoints=2., borderpad=0.,
                        handletextpad=0.2, loc=2)

    for handle in lgd.legendHandles:
        handle._legmarker.set_markersize(3.)
        handle.set_linewidth(1.)

    # save the figure
    if var == 'Vcmax':
        namefig = os.path.join(get_fig_dir(), 'kmax_2_Vcmax')

    elif var == 'Tair':
        namefig = os.path.join(get_fig_dir(), 'kmax_2_Tair')

    elif var == 'both':
        plt.tight_layout()
        figdir = os.path.join(get_fig_dir(), 'final_4_paper')
        namefig = os.path.join(figdir, 'profitmax_2_forcing')

    fig.savefig('%s.png' % (namefig), dpi=1000, bbox_inches='tight')
    fig.savefig('%s.pdf' % (namefig), dpi=1800, bbox_inches='tight')
    plt.close()

    return


#=======================================================================

if __name__ == "__main__":

    # run mode
    res = 'high'
    photo = 'Farquhar'
    case = 2
    parallel = False
    solstep = 'var'
    symbolic = False

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.weight'] = 'light'
    plt.rcParams['axes.titlepad'] = 30

    # purple, blue, green, orange, gray
    colours = ['#762a83', '#4393c3', '#7fbc41', '#fdb462', 'gray']

    try:
        plt.rc('text', usetex=True)

        try:
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
                                                   r'\usepackage{cmbright}']
            main()

        except Exception as e:
            try:
                plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
                main()

            except Exception as e:
                plt.rcParams['text.latex.preamble'] = \
                    [r'\usepackage{cmbright}']
                main()

    except Exception as e:
        main()
