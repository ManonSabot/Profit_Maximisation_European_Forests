#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions used to calculate and plot the effective g1 values predicted
by the Profitmax model and their sensitivities to VPD

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "plot the site-level effective g1 and VPD sensitivities"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (26.10.2019)"
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
from scipy.optimize import curve_fit  # fit models to one another

# plotting modules
import matplotlib.pyplot as plt
import string  # automate subplot lettering
from matplotlib.ticker import MaxNLocator
import matplotlib.lines as mlines  # for the legend

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM import conv, cst  # unit converter & general constants
from TractLSM.Utils import read_csv  # read in data files
from TractLSM import dparams  # model parameters
from TractLSM.run_utils import time_step  # read each time step's params
from TractLSM.SPAC import net_radiation, leaf_energy_balance, vpsat
from TractLSM.SPAC import hydraulics, leaf_temperature, calc_photosynthesis
from TractLSM.CH2OCoupler.ProfitMax import hydraulic_cost, photo_gain
from TractLSM.CH2OCoupler.USO import calc_trans  # needed USO model

from plot_utils import get_best_kmax_calib  # best kmax's number name


#=======================================================================

def main():

    """
    Main: Explores the relative behaviours of the effective g1 given by
          the ProfitMax compared to those from the USO model, as well as
          their varying sensitivities to VPD.

    Returns:
    --------
    output/figures/final_4_paper/sensitivities_to_VPD.png

    """

    # read the site params
    main_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

    while 'src' in main_dir:
        main_dir = os.path.dirname(main_dir)

    fluxsites_folder = os.path.join(os.path.join(main_dir, 'input'),
                                    'fluxsites')

    # create artificial forcing data
    df = make_artificial_data()

    # get the site info
    info, __ = read_csv(os.path.join(fluxsites_folder, 'info_sites.csv'))
    info.index = info['Site']

    # run the big-leaf ProfitMax model and fit a g1 and sensitivity to D
    g1, s = Sperry_2_Medlyn(df, info)

    # colours (purple, blue, green, yellow, orange, pink)
    col = ['#762a83', '#bcbddc', '#2c7fb8', '#014636', '#7fbc41', '#dfc27d',
           '#7f2704', '#fdae61', '#e31a1c', '#de77ae']

    plt.figure(figsize=(8., 3.5))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    plt.subplots_adjust(wspace=0.125)

    # reference Medlyn g1
    sites = sorted(list(np.unique(np.array(info.index.values))))
    g1_ref = [info.loc[site, 'g1'] for site in sites]

    # plot the g1 data
    ax1.scatter(g1 ** (0.5 / s), s, c=col)
    ax1.scatter(g1_ref, s, marker='^', facecolor='none', edgecolor=col)

    for i in range(len(s)):  # lines between the two g1 values

        ax1.plot([np.minimum(g1[i] ** (0.5 / s[i]), g1_ref[i]),
                 np.maximum(g1[i] ** (0.5 / s[i]), g1_ref[i])],
                 [s[i], s[i]], linestyle=':', color=col[i])

    print('average sensitivity : %s, median sensitivity : %s'
          % (str(round(np.mean(s), 2)), str(round(np.median(s), 2))))

    # effects of the different sensitivities, for a given g1
    df, sensitivities = gs_sensitivities(sites, s)
    tr = 0

    for which in sensitivities:

        x = np.log(df.loc[which, 'D_all'])
        y = df.loc[which, 'gs_all']

        if which == 'Ref.':  # sensitivity of 0.5
            ax2.plot(x, y, ':k', linewidth=4., label=which, zorder=20)

        else:  # site-level sensitivities
            ax2.plot(x, y, col[tr], label=which)
            tr += 1

    # restrict the number of ticks
    ax1.xaxis.set_major_locator(MaxNLocator(5))
    ax1.yaxis.set_major_locator(MaxNLocator(5))
    ax2.set_xticks([np.log(0.1), np.log(0.2), np.log(0.4), np.log(0.8),
                   np.log(1.6), np.log(3.2)])
    ax2.set_xticklabels(['0.1', '0.2', '0.4', '0.8', '1.6', '3.2'])
    ax2.yaxis.set_major_locator(MaxNLocator(5))

    # number the subplots
    ax1.text(0.92, 0.95, '(%s)' % (string.ascii_lowercase[0]),
             transform=ax1.transAxes, fontweight='heavy', fontsize=10)
    ax2.text(0.92, 0.95, '(%s)' % (string.ascii_lowercase[1]),
             transform=ax2.transAxes, fontweight='heavy', fontsize=10)

    # axes labels
    ax1.set_xlabel('g$_{1}$ (kPa$^{0.5}$)')
    ax1.set_ylabel('$\sigma$')
    ax2.set_xlabel('D (kPa)')
    ax2.set_ylabel('g$_{s}$ (mol m$^{-2}$ s$^{-1}$)')

    # add legend
    handles, labels = ax2.get_legend_handles_labels()
    handles.append(handles.pop(labels.index('Ref.')))
    labels.append(labels.pop(labels.index('Ref.')))

    # add the marker symbols at beginning of legend
    circle = mlines.Line2D([], [], marker='o', color='k', linestyle='none')
    triangle = mlines.Line2D([], [], marker='^', markeredgecolor='k',
                             markerfacecolor='none', linestyle='none')
    handles = [circle, triangle] + handles
    labels = [r'Profit$_{\rm max}$ g$_{1}$', 'Control g$_{1}$'] + labels
    plt.tight_layout()
    ax2.legend(handles, labels, bbox_to_anchor=(1.025, 0.98))

    # save the figure
    figdir = os.path.join(get_fig_dir(), 'final_4_paper')
    plt.savefig(os.path.join(figdir, 'sensitivities_to_VPD.png'), dpi=1000,
                bbox_inches='tight')
    plt.savefig(os.path.join(figdir, 'sensitivities_to_VPD.eps'), dpi=600,
                bbox_inches='tight')
    plt.close()

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
    p = p.append(pd.Series([12., 25., 1., 2000., 1., 1., p.albedo_ws],
                 index=['hod', 'Tair', 'VPD', 'PPFD', 'coszen',
                        'scale2can', 'albedo_s']))

    return p


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


def make_artificial_data():

    """
    Generates a dataframe of artificial data, with co-varying VPD and
    air temperature (bound by RH between 45 and 95%).

    Returns:
    --------
    df: pandas dataframe
        dataframe containing the forcing data

    """

    p = declared_params()
    VPD = np.arange(0.05, 3.05, 0.05)
    df = pd.DataFrame(np.tile(p, (len(VPD), 1)), index=range(0, len(VPD)),
                      columns=p.index)
    df['VPD'] = VPD
    df['hod'] = np.tile([12.5, 12.], int(len(VPD) / 2.))

    # bound the temperature for it to be valid
    RH = np.arange(45., 95., 5.)
    df = pd.DataFrame(np.tile(df, (len(RH), 1)),
                      index=range(0, len(VPD) * len(RH)), columns=df.columns)
    df['RH'] = 1.

    for i in range(len(RH)):

        df.loc[len(VPD)*i:len(VPD)*(i+1)-1, 'Tair'] = Tair_from_VPD(VPD, RH[i])
        df.loc[len(VPD)*i:len(VPD)*(i+1)-1, 'RH'] = RH[i]

    df.fillna(0., inplace=True)

    return df


def kmax_values(sites):

    """
    Retrieves the site values of the 'best' kmax from the calibration
    project.

    Arguments:
    ----------
    sites: array
        names of the sites to plot

    Returns:
    --------
    kmax: array
        maximum hydraulic conductance [mmol m-2 s-1 MPa-1]

    """

    basedir = os.path.dirname(os.path.realpath(sys.argv[0]))

    while 'src' in basedir:
        basedir = os.path.dirname(basedir)

    iname = os.path.join(os.path.join(basedir, 'input'), 'projects')
    iname = os.path.join(iname, 'var_kmax_sample')

    basedir = os.path.join(os.path.join(os.path.join(basedir, 'output'),
                           'projects'), 'var_kmax')
    kmax = get_best_kmax_calib(sites, basedir)
    sites = [s + kmax[i] for i, s in enumerate(sites)]

    kmax = np.zeros(len(sites))

    for i in range(len(sites)):

        for file in os.listdir(iname):

            if ((file.endswith('.csv')) and ('actual' in file) and
               ('%s_' % (sites[i]) in file)):
                df, __ = read_csv(os.path.join(iname, file))
                kmax[i] = df['kmax'].values[0]

    return kmax


def update_site_params(df, p):

    """
    Updates a df with site-specific parameters.

    Arguments:
    ----------
    df: pandas dataframe
        dataframe containing all the data

    p: pandas series
        site parameters

    Returns:
    --------
    df: pandas dataframe
        updated model's parameters (with site info)

    """

    try:
        if str(p.Vmax25) != str(pd.np.NaN):
            df['Vmax25'] = p.Vmax25

    except AttributeError:
        pass

    try:
        if str(p.albedo_l) != str(pd.np.NaN):
            df['albedo_l'] = p.albedo_l

    except AttributeError:
        pass

    try:
        if str(p.max_leaf_width) != str(pd.np.NaN):
            df['max_leaf_width'] = p.max_leaf_width

    except AttributeError:
        pass

    try:
        if str(p.P50) != str(pd.np.NaN):
            df['P50'] = p.P50

    except AttributeError:
        pass

    try:
        if str(p.P88) != str(pd.np.NaN):
            df['P88'] = p.P88

    except AttributeError:
        pass

    try:
        if str(p.ratiocrit) != str(pd.np.NaN):
            df['ratiocrit'] = p.ratiocrit

    except AttributeError:
        pass

    try:
        if str(p.Psie) != str(pd.np.NaN):
            df['Psie'] = p.Psie
            df['Ps'] = p.Psie

    except AttributeError:
        pass

    try:
        if str(p.g1) != str(pd.np.NaN):
            df['g1'] = p.g1

    except AttributeError:
        pass

    try:
        if str(p['Composite LAI']) != str(pd.np.NaN):
            df['LAI'] = p['Composite LAI']

    except AttributeError:
        pass

    return df


def profit_psi(p, photo='Farquhar', res='low', case=2):

    """
    Simplified big leaf that finds the instateneous profit maximization,
    following the optimization criterion for which, at each instant in
    time, the stomata regulate canopy gas exchange and pressure to
    achieve the maximum profit, which is the maximum difference between
    the normalized photosynthetic gain (gain) and the hydraulic cost
    function (cost). That is when d(gain)/dP = d(cost)/dP.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    res: string
        either 'low' (default), 'med', or 'high' to run the optimising
        solver

    case: int
        either 1, where Ci isn't stricktly bound by the standard model,
        or 2 where it is

    Returns:
    --------
    gs: float
        stomatal conductance [mol m-2 s-1] at maximum profit

    An: float
        net photosynthetic assimilation rate [umol m-2 s-1] at maximum
        profit

    Dleaf: float
        leaf vapour pressure deficit [kPa] at maximum profit

    """

    # net radiation
    p.Rnet = net_radiation(p)

    esat_a = vpsat(p.Tair)  # kPa

    # big leaf scaling (eq 4; Haxeltine & Prentice 1996)
    scale_An = 1. - np.exp(-0.5 * p.LAI)  # applied on APAR
    scale_E = 2. * (1. - np.exp(-0.5 * p.LAI))  # integrated over APAR

    # hydraulics
    P, trans = hydraulics(p, res=res)
    cost, __ = hydraulic_cost(p, P)

    # apply the scaling on E
    trans *= scale_E

    try:
        gain, Ci = photo_gain(p, trans, photo, res, case, parallel, solstep,
                              symbolic)

        # look for the most net profit
        profit = gain - cost

        # deal with edge cases by rebounding the solution
        gc, gs, __, __ = leaf_energy_balance(p, trans)
        profit_check = profit[1:][(gc[1:] > cst.zero)]
        idx = np.isclose(profit, max(profit_check))
        idx = [list(idx).index(e) for e in idx if e]

        if idx:  # opt values
            gs = gs[idx[0]]

            # rubisco- or electron transport-limitation?
            p.PPFD *= scale_An
            An, __, __ = calc_photosynthesis(p, trans[idx[0]], Ci[idx[0]],
                                             photo)
            p.PPFD /= scale_An

            # saturation vapour pressure of water at T
            Tleaf, __ = leaf_temperature(p, trans[idx[0]])
            esat_l = vpsat(Tleaf)  # kPa
            Dleaf = (esat_l - (esat_a - p.VPD))  # leaf-air vpd, kPa

            # if the critical point has been reached, stall
            if np.isclose(trans[idx[0]], 0.):
                gs, An, Dleaf = (0.,) * 3

            if any(np.isnan([trans[idx[0]], gs, An, Dleaf])):
                gs, An, Dleaf = (0.,) * 3

    except (ValueError, TypeError):  # no opt
        gs, An, Dleaf = (0.,) * 3

    return gs, An, Dleaf


def solve_std(p, photo='Farquhar', s=0.5, threshold_conv=0.015, iter_max=40):

    """
    Simplified big leaf that checks the energy balance by looking for
    convergence of the new leaf temperature with the leaf temperature
    predicted by the previous iteration. Then returns the corresponding
    gs, An, Dleaf.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    threshold_conv: float
        convergence threshold for the new leaf temperature to be in
        energy balance

    iter_max: int
        maximum number of iterations allowed on the leaf temperature
        before reaching the conclusion that the system is not energy
        balanced

    Returns:
    --------
    gs: float
        stomatal conductance [mol m-2 s-1]

    An: float
        net photosynthetic assimilation rate [umol m-2 s-1]

    Dleaf: float
        leaf vapour pressure deficit [kPa]

    """

    # net radiation
    p.Rnet = net_radiation(p)

    esat_a = vpsat(p.Tair)  # kPa

    # big leaf scaling (eq 4; Haxeltine & Prentice 1996)
    scale_An = 1. - np.exp(-0.5 * p.LAI)  # applied on APAR
    scale_E = 2. * (1. - np.exp(-0.5 * p.LAI))  # integrated over APAR

    # initialisation
    Cs = p.CO2  # Pa
    Tleaf = p.Tair  # deg C
    Dleaf = np.maximum(0.05, p.VPD)  # gs model not valid < 0.05

    # initialise gs_over_A
    g0 = 1.e-9  # g0 ~ 0, removing it entirely introduces errors
    Cs_umol_mol = Cs * conv.MILI * conv.FROM_kPa  # umol mol-1
    gs_over_A = g0 + (1. + p.g1 / (Dleaf ** s)) / Cs_umol_mol

    # iter on the solution until it is stable enough
    iter = 0

    while True:

        p.PPFD *= scale_An
        An, __, __, __ = calc_photosynthesis(p, 0., Cs, photo, Tleaf=Tleaf,
                                             gs_over_A=gs_over_A)
        p.PPFD /= scale_An

        # stomatal conductance
        Cs_umol_mol = Cs * conv.MILI * conv.FROM_kPa
        gs_over_A = (g0 + (1. + p.g1 / (Dleaf ** s)) / Cs_umol_mol)
        gs = np.maximum(cst.zero, conv.GwvGc * gs_over_A * An)

        # calculate new trans, gw, gb, Tleaf
        trans, real_zero, gw, gb = calc_trans(p, Tleaf, gs)
        trans *= scale_E
        new_Tleaf, __ = leaf_temperature(p, trans, Tleaf=Tleaf)

        # new Cs (in Pa)
        boundary_CO2 = (conv.ref_kPa * conv.FROM_MILI * An /
                        (gb * conv.GbcvGb + gs * conv.GcvGw))
        Cs = np.maximum(cst.zero, np.minimum(p.CO2, p.CO2 - boundary_CO2))

        if (np.isclose(trans, cst.zero, rtol=cst.zero, atol=cst.zero) or
           np.isclose(gw, cst.zero, rtol=cst.zero, atol=cst.zero) or
           np.isclose(gs, cst.zero, rtol=cst.zero, atol=cst.zero)):
            Dleaf = np.maximum(0.05, p.VPD)  # kPa

        else:
            esat_l = vpsat(new_Tleaf)  # kPa
            Dleaf = (esat_l - (esat_a - p.VPD))  # leaf-air vpd, kPa

        # force stop when atm. conditions yield E < 0. (non-physical)
        if (iter < 1) and (not real_zero):
            real_zero = None

        # check for convergence
        if ((real_zero is None) or (iter > iter_max) or ((real_zero) and
           (abs(Tleaf - new_Tleaf) <= threshold_conv) and not
           np.isclose(gs, cst.zero, rtol=cst.zero, atol=cst.zero))):
            break

        # no convergence, iterate on leaf temperature
        Tleaf = new_Tleaf
        iter += 1

    # if the critical point has been reached, stall
    if np.isclose(trans, 0.):
        gs, An = (0.,) * 2

    if any(np.isnan([trans, gs, An, Dleaf])):
        gs, An, Dleaf = (0.,) * 3

    return gs, An, Dleaf


def fit_gs(X, a, b):

    """
    Function to fit gs to: this is used to find the parameters a and b.

    Arguments:
    ----------
    X: tuple
        contains the forcing variables x (i.e. VPD, [kPa]) and y
        (i.e. A/Ca, [mol m-2 s-1])

    a: float
        g1 parameter which must be fitted [kPa^sensitivity]

    b: float
        sensitivity parameter which must be fitted [unitless]

    Returns:
    --------
    gs: float
        stomatal conductance [mol m-2 s-1]

    """

    x, y = X

    return 1.e-9 + conv.GwvGc * (1. + a * x ** (-b)) * y


def fit_2_Medlyn(gs, D, AoCa, df):

    """
    Fits predictions of stomatal conductance to a model based on VPD and
    the ratio of AoCa.

    Arguments:
    ----------
    gs: array
        predicted stomatal conductance [mol m-2 s-1]

    D: array
        predicted leaf level vapour pressure deficit [kPa]

    AoCa: array
        predicted ratio of An to Ca [mol m-2 s-1]

    df: pandas dataframe
        dataframe containing the forcing data

    Returns:
    --------
    g1: float
        plant water use [kPa^s]

    s: float
        stomatal senstivity to VPD [unitless]

    """

    # find g1: @ D = 1 kPa, g1 = gs * Ca / (A * 1.57) - 1
    mask = np.logical_and(np.isclose(D, 1., rtol=5.e-02, atol=5.e-02),
                          np.logical_and(df.Tair > 20., df.Tair < 30.))
    g1 = np.median(gs[mask] / (AoCa[mask] * conv.GwvGc) - 1.)

    # weights for the fits to be sensible (because of the RH variations)
    Dbins = np.digitize(D, np.unique(df['VPD']))
    __, weights = np.unique(Dbins, return_counts=True)
    idx = np.cumsum(weights)
    weights = weights.astype(float) / float(np.amax(weights))
    weights = ([np.asarray([weights[0]])] +
               [np.ones(len(df))[idx[i]:idx[i+1]] * weights[i] for
                i in range(len(idx) - 1)] + [np.asarray([weights[-1]])])
    weights = np.concatenate(weights).ravel()

    # make sure we have the right number of weights
    while len(weights) > len(D):
        weights = weights[1:]

    while len(weights) < len(D):
        weights = np.insert(weights, 0, weights[0])

    weights /= np.sum(weights)  # normalised total to 1.

    # bins to fit
    Dbins = ([D[Dbins == i] for i in range(1, len(np.unique(df['VPD'])) + 1)]
             [1:])

    # now, optimise the curve fitting
    s = 0.
    total_weight = 0.
    p0 = 0.5  # initial guess for the gs-D sensitivity

    for abin in Dbins:

        try:
            if len(abin) > 3:
                mask = np.logical_and(D >= np.amin(abin), D <= np.amax(abin))

                X = D[mask], AoCa[mask]
                params, __ = curve_fit(lambda x, b: fit_gs(X, g1, b), X,
                                       gs[mask], p0, method='trf')

                if params[0] > 0.:
                    s += params[0] * np.sum(weights[mask])
                    total_weight += np.sum(weights[mask])

        except (RuntimeError, TypeError):
            pass

    s /= total_weight

    return g1, s


def Sperry_2_Medlyn(df, info):

    """
    Fits the gs predicted by the ProfitMax algorithm to the gs predicted
    by the USO algorithm across sites, by adjusting the g1 and
    sensitivity terms in USO.

    Arguments:
    ----------
    df: pandas dataframe
        dataframe containing the forcing data

    info: pandas dataframe
        dataframe containing all the site specific parameters

    Returns:
    --------
    g1: array
        plant water use [kPa^s]

    s: array
        stomatal senstivity to VPD [unitless]

    """

    # get the site calibrated kmax
    sites = sorted(list(np.unique(np.array(info.index.values))))
    kmax = kmax_values(sites)

    # declare empty arrays to store the g1 and s values
    g1 = np.zeros(len(sites))
    s = np.zeros(len(sites))

    tr = 0

    for site in sites:

        # restrict on per site basis to avoid cold T effects in An model
        Tmin = np.minimum(10., info.loc[site, 'CRU cTair'] - 5.)
        Tmax = np.maximum(40., info.loc[site, 'CRU cTxx'] + 15.)
        mask = np.logical_and(df['Tair'] > Tmin, df['Tair'] < Tmax)
        df1 = df.copy()[mask]
        df1.reset_index(inplace=True, drop=True)
        df1.sort_values(by=['VPD'], inplace=True)
        df1.reset_index(inplace=True, drop=True)

        # update site parameters
        df1 = update_site_params(df1, info.loc[site])
        df1.kmax = kmax[tr]

        gs_psi = np.zeros(len(df1))
        A_psi = np.zeros(len(df1))
        D_psi = np.zeros(len(df1))

        for i in range(len(df1)):

            p = time_step(df1, i)
            gs_psi[i], A_psi[i], D_psi[i] = profit_psi(p, photo=photo, res=res)

        AoCa = A_psi / (p.CO2 * 1.e3 / 101.325)
        g1[tr], s[tr] = fit_2_Medlyn(gs_psi, D_psi, AoCa, df1)

        print('%s g1: %s, sensitivity: %s' % (site, str(round(g1[tr], 2)),
              str(round(s[tr], 2))))

        tr += 1

    return g1, s


def gs_sensitivities(sites, sens):

    """
    Calculates the gs predicted by the USO algorithm for a g1 varying
    with the site-level sensitivities to VPD.

    Arguments:
    ----------
    sites: array
        names of the sites to plot

    sens: array
        stomatal senstivity to VPD [unitless]

    Returns:
    --------
    df2: pandas dataframe
        dataframe containing the predicted gs for each value of sens

    sens: array
        stomatal senstivity to VPD [unitless] in the order that matches
        df2

    """

    # generate very simple forcing data, where VPD only varies
    p = declared_params()
    VPD = np.arange(0.05, 3.05, 0.010)
    df1 = pd.DataFrame(np.tile(p, (len(VPD), 1)), index=range(0, len(VPD)),
                       columns=p.index)
    df1['VPD'] = VPD
    df1['hod'] = np.tile([12.5, 12.], int(len(VPD) / 2.))
    df1['Tair'] = 25.  # impose 25 deg C at anytime
    df1.drop_duplicates(inplace=True)
    df1.reset_index(inplace=True, drop=True)
    df1.sort_values(by=['VPD'], inplace=True)
    df1.reset_index(inplace=True, drop=True)

    g1 = 2.  # make g1 = 2 @ sens = 0.5
    df1['g1'] = g1

    df2 = pd.DataFrame(np.nan, index=np.arange((len(sites) + 1) * len(df1)),
                       columns=['VPD', 'sensitivity', 'gs_all', 'D_all'])

    # run the USO model with the standard sensitivity
    for i in range(len(df1)):

        p = time_step(df1, i)
        df2.loc[i, 'VPD'] = p.VPD
        df2.loc[i, 'sensitivity'] = 'Ref.'
        df2.loc[i, 'gs_all'], A, df2.loc[i, 'D_all'] = solve_std(p,
                                                                 photo=photo)

    tr = 1

    # run the USO model with the site-level sensitivities
    for site in sites:

        df1['g1'] = g1 ** (0.5 / sens[tr - 1])  # adjust g1

        for i in range(len(df1)):

            p = time_step(df1, i)
            j = tr * len(df1) + i
            df2.loc[j, 'VPD'] = p.VPD
            df2.loc[j, 'sensitivity'] = '%s' % (site)
            df2.loc[j, 'gs_all'], A, df2.loc[j, 'D_all'] = \
                solve_std(p, photo=photo, s=sens[tr - 1])

        tr += 1

    sens = np.unique(df2['sensitivity'])

    # average the fits per VPD value
    df2 = df2.groupby(['sensitivity', 'VPD']).max()
    df2.index.droplevel(level=1)

    return df2, sens


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


#=======================================================================

if __name__ == "__main__":

    # run mode
    res = 'low'
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
