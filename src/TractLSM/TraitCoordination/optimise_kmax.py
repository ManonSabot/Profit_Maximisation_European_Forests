# -*- coding: utf-8 -*-

"""
Hypothetising long-term coordination between the photosynthetic traits,
the hydraulic traits and the climate to calculate maximum hydraulic
conductance.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Haxeltine, Alex, and I. Colin Prentice. "BIOME3: An equilibrium
  terrestrial biosphere model based on ecophysiological constraints,
  resource availability, and competition among plant functional types."
  Global Biogeochemical Cycles 10.4 (1996): 693-709.

"""

__title__ = "Optimise maximum hydraulic conductance"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (20.08.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import numpy as np  # array manipulations, math operators
import bottleneck as bn  # faster C-compiled np for all nan operations
import pandas as pd  # read/write dataframes, csv files
from scipy.integrate import quad  # integrate on a range
from sympy import Symbol  # express functions symbolically
from sympy.utilities.lambdify import lambdify  # expression to function

# own modules
from TractLSM import conv, cst  # unit converter & general constants
from TractLSM.SPAC import f, Weibull_params  # hydraulics
from TractLSM.SPAC import calc_photosynthesis  # photosynthesis
from TractLSM.CH2OCoupler.ProfitMax import A_trans, hydraulic_cost


# ======================================================================

class df2class:

    def __init__(self, df):

        """
        Converts a dataframe or dataseries to a class.

        """

        for name, value in df.items():
            setattr(self, name, value)

        return


def kmax_matrices(opt1, opt2, opt3, fun, p, strategy, P50, P88, ratiocrit):

    """
    Calculates three possible values of maximum hydraulic conductance,
    each of which corresponds to a different optimal point on the
    hydraulic stream.

    Arguments:
    ----------
    opt1: array
        empty data array

    opt2: array
        empty data array

    opt3: array
        empty data array

    fun: function
        function of Ci, the integrated vulnerability curve, and the
        optimal sol

    p: pandas series
        model's parameters

    strategy: string
        'all', 'high' (at P12), 'optimal', or 'low' (before the cost
        offsets the profit) solution

    P50: array
        leaf water potential [MPa] causing a 50% decline in hydraulic
        conductance

    P88: array
        leaf water potential [MPa] causing a 88% decline in hydraulic
        conductance

    ratiocrit: array
        ratio setting the critical lowest possible hydraulic
        conductance, indicative of how conservative is the stomatal
        behaviour

    Returns:
    --------
    opt1: array
        data array containing solution1 and the specific P50, P88,
        rcrit, Vmax25, VPD, and CO2 values used to solve

    opt2: array
        data array containing solution2 and the specific P50, P88,
        rcrit, Vmax25, VPD, and CO2 values used to solve

    opt3: array
        data array containing solution2 and the specific P50, P88,
        rcrit, Vmax25, VPD, and CO2 values used to solve

    """

    for i, j in zip(range(len(P50)), range(len(P88))):  # combis
        p.P50 = P50[i]
        p.P88 = P88[j]

        if p.P88 > p.P50:
            b, c = Weibull_params(p)  # MPa, unitless

            # the optimal kmax can only be calculated if c >= 0.9
            if c < 0.9:  # exponential (open vessel artifact?)
                c = 1.  # assume 1., safer than 0.9
                b = p.P50 / ((- np.log(1 - 0.5)) ** (1. / c))

            for k in range(len(ratiocrit)):
                rcrit = ratiocrit[k]

                # get the downstream leaf water potentials P (MPa)
                Pcrit = - b * np.log(1. / rcrit) ** (1. / c)  # MPa
                P = np.linspace(p.Ps, Pcrit, 2000)  # from Ps to Pcrit
                int_vuln = np.empty_like(P)

                # get P12
                P12 = - b * ((- np.log(1 - 0.12)) ** (1. / c))

                for l in range(len(P)):

                    int_vuln[l], err = quad(f, P[l], P[0], args=(b, c))

                kmax = np.linspace(rcrit, 25., 20000)
                sols = fun(0.7 * p.CO2, np.expand_dims(int_vuln, axis=1), kmax)
                idx = bn.nanargmin(abs(sols), axis=1)  # sols ~0.

                # solution with approximate minimizing
                sols = np.asarray([kmax[idx[e]] for e in range(len(P))])[1:]

                if (strategy == 'all') or (strategy == 'high'):
                    sol1 = cst.zero

                    try:  # P12 sol
                        sol1 = sols[bn.nanargmin(P[1:] - P12 >= 0.) - 1]

                        if sol1 == max(kmax):
                            sol1 = cst.zero

                        opt1[i, j, k, :] = [sol1, p.P50, p.P88, rcrit,
                                            p.Vmax25, p.VPD, p.CO2]

                    except ValueError:
                        pass

                if ((strategy == 'all') or (strategy == 'optimal') or
                   (strategy == 'low')):
                    profit_max = -1.  # initialise opt maximal profit
                    gain_max = 0.  # initialise low maximal gain
                    sol2 = cst.zero
                    sol3 = cst.zero

                    for temp_sol in sols:  # loop over possible opt kmax

                        trans = temp_sol * conv.FROM_MILI * int_vuln[1:]
                        A_P = A_trans(p, trans, 0.7 * p.CO2, Tleaf=p.Tref)
                        gain = A_P / np.amax(A_P)

                        if np.amax(A_P) < 0.:  # resp >> An everywhere
                            gain *= -1.

                        cost, __ = hydraulic_cost(p, P[1:])
                        profit = (gain - cost)[:-1]

                        if (strategy == 'all') or (strategy == 'optimal'):

                            try:
                                if profit[bn.nanargmax(profit)] > profit_max:
                                    sol2 = sols[bn.nanargmax(profit)]
                                    profit_max = profit[bn.nanargmax(profit)]

                            except ValueError:
                                sol2 = cst.zero

                            if sol2 == max(kmax):
                                sol2 = cst.zero

                            opt2[i, j, k, :] = [sol2, p.P50, p.P88, rcrit,
                                                p.Vmax25, p.VPD, p.CO2]

                        if (strategy == 'all') or (strategy == 'low'):

                            try:
                                ilow = (bn.nanargmin(profit - cost[:-1] > 0.) -
                                        1)

                                try:
                                    if profit[ilow] > gain_max:
                                        sol3 = sols[ilow]
                                        gain_max = profit[ilow]

                                except ValueError:
                                    pass

                            except ValueError:
                                pass

                            if sol3 == max(kmax):
                                sol3 = cst.zero

                            opt3[i, j, k, :] = [sol3, p.P50, p.P88, rcrit,
                                                p.Vmax25, p.VPD, p.CO2]

    return opt1, opt2, opt3


def optimal_kmax(p, photo, strategy='all', P50=None, P88=None, ratiocrit=None,
                 Vcmax25=None, VPD=None, CO2=None):

    """
    Calculates the assimilation rate given the supply function, gc. No
    respiration here as this is the "physics" An. Both An and E are
    scaled up following Haxeltine & Prentice (1996).

    Arguments:
    ----------
    p: pandas series
        model's parameters

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    strategy: string
        'all', 'high' (at P12), 'optimal', or 'low' (before the cost
        offsets the profit) solution

    P50: float or array
        leaf water potential [MPa] causing a 50% decline in hydraulic
        conductance

    P88: float or array
        leaf water potential [MPa] causing a 88% decline in hydraulic
        conductance

    ratiocrit: float or array
        ratio setting the critical lowest possible hydraulic
        conductance, indicative of how conservative is the stomatal
        behaviour

    Vcmax25: float or array
        maximum carboxylation rate at 25 degC [umol.m-2.s-1]

    VPD: float or array
        vapour pressure deficit [kPa]

    CO2: float or array
        atmospheric CO2 concentration [Pa]

    Returns:
    --------
    Three data arrays, each containing a solution for kmax and the
    specific P50, P88, rcrit, Vmax25, VPD, and CO2 values used to solve.

    """

    # archive the parameter class before it gets tweaked
    archive = p.copy()

    # pass the df to class to allow for symbols in it!
    p = df2class(p)

    # big leaf scaling (eq 4; Haxeltine & Prentice 1996)
    scale_An = 1. - np.exp(-0.5 * p.LAI)  # applied on APAR
    scale_E = 2. * (1. - np.exp(-0.5 * p.LAI))  # integrated over APAR

    # declare symbols here
    CI = Symbol('CI', positive=True)
    FP = Symbol('FP', positive=True)  # integrated vulnerability curve
    kmax = Symbol('kmax', positive=True)
    p.Vmax25 = Symbol('Vcmax', positive=True)
    p.CO2 = Symbol('CO2', positive=True)
    p.VPD = Symbol('VPD', positive=True)
    E = kmax * conv.FROM_MILI * FP * scale_E
    p.Rlref = 0.015 * p.Vmax25

    # get expression of the optimal kmax
    A_P = A_trans(p, E, CI, Tleaf=p.Tair)
    p.PPFD *= scale_An
    A_Ci, __, __ = calc_photosynthesis(p, E, CI, photo, Tleaf=p.Tair)

    expr = A_Ci - A_P

    if P50 is None:
        P50 = archive.P50

    try:
        len(P50)

    except TypeError:
        P50 = [P50]

    if P88 is None:
        P88 = archive.P88

    try:
        len(P88)

    except TypeError:
        P88 = [P88]

    if ratiocrit is None:
        ratiocrit = archive.ratiocrit

    try:
        len(ratiocrit)

    except TypeError:
        ratiocrit = [ratiocrit]

    if Vcmax25 is None:
        Vcmax25 = archive.Vmax25

    try:
        len(Vcmax25)

    except TypeError:
        Vcmax25 = [Vcmax25]

    if VPD is None:
        VPD = archive.VPD

    try:
        len(VPD)

    except TypeError:
        VPD = [VPD]

    if CO2 is None:
        CO2 = archive.CO2

    try:
        len(CO2)

    except TypeError:
        CO2 = [CO2]

    # declare size of solution matrices
    Kmax_opt1 = np.ones((len(P50), len(P88), len(ratiocrit), len(Vcmax25),
                        len(VPD), len(CO2), 7))
    Kmax_opt2 = np.empty_like(Kmax_opt1)
    Kmax_opt3 = np.empty_like(Kmax_opt1)

    # make the parameter class a pandas series again
    attrs = vars(p)
    p = {item[0]: item[1] for item in attrs.items()}
    p = pd.Series(p)

    # update the soil water potential with the input info
    if archive.Ps > archive.Psie:
        p.Ps = p.Psie

    else:
        p.Ps = archive.Ps

    for dim1 in range(len(Vcmax25)):

        for dim2 in range(len(VPD)):

            for dim3 in range(len(CO2)):

                # redeclare symbols to allow solving
                p.Vmax25 = Symbol('Vcmax', positive=True)
                p.CO2 = Symbol('CO2', positive=True)
                p.VPD = Symbol('VPD', positive=True)

                # generate function ~ 0. for opt kmax from expression
                new_expr = expr.subs({p.Vmax25: Vcmax25[dim1],
                                      p.VPD: VPD[dim2],
                                      p.CO2: CO2[dim3]})

                fun = lambdify((CI, FP, kmax), new_expr, 'numpy')

                # update Vcmax, VPD and CO2
                p.Vmax25 = Vcmax25[dim1]
                p.VPD = VPD[dim2]
                p.CO2 = CO2[dim3]
                p.Rlref = 0.015 * p.Vmax25

                # fill the matrices with all the potential combinations
                (Kmax_opt1[:, :, :, dim1, dim2, dim3, :],
                 Kmax_opt2[:, :, :, dim1, dim2, dim3, :],
                 Kmax_opt3[:, :, :, dim1, dim2, dim3, :]) = \
                    (kmax_matrices(Kmax_opt1[:, :, :, dim1, dim2, dim3, :],
                                   Kmax_opt2[:, :, :, dim1, dim2, dim3, :],
                                   Kmax_opt3[:, :, :, dim1, dim2, dim3, :],
                                   fun, p, strategy, P50, P88, ratiocrit))

    return np.squeeze(Kmax_opt1), np.squeeze(Kmax_opt2), np.squeeze(Kmax_opt3)
