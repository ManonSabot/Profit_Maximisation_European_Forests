# -*- coding: utf-8 -*-

"""
Functions related to plant hydraulics, used to calculate xylem
vulnerability to cavitation, hydraulic conductance, and transpiration.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Neufeld et al. (1992). Genotypic variability in vulnerability of leaf
  xylem to cavitation in water-stressed and well-irrigated sugarcane.
  Plant Physiology, 100(2), 1020-1028.
* Sperry et al. (2017). Predicting stomatal responses to the environment
  from the optimization of photosynthetic gain and hydraulic cost.
  Plant, cell & environment, 40(6), 816-830.

"""

__title__ = "Plant hydraulics module"
__author__ = "Manon E. B. Sabot"
__version__ = "2.0 (29.11.2017)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import numpy as np  # array manipulations, math operators
from scipy.integrate import quad  # integrate on a range

# own modules
from TractLSM import conv, cst  # unit converter & general constants


# ======================================================================

def f(P, b, c):

    """
    Calculates soil and xylem vulnerability curves based on a
    two-parameter Weibull function (Neufeld et al., 1992). In turn, this
    is used to calculate hydraulic conductance and transpiration at
    steady-state.

    Arguments:
    ----------
    P: array
        leaf or soil (upstream) water potential [MPa], an array of
        values from the soil water potential Ps to the critical water
        potential Pcrit for which cavitation of the xylem occurs

    b: float
        one of two Weibull parameters [MPa], it is P at k / kmax = 0.37

    c: float
        one of two Weibull parameters [unitless], it controls whether
        the curve has a threshold sigmoidal form or non-threshold
        sigmoidal form

    Returns:
    --------
    The vulnerability curves of the the plant [unitless] describing the
    response of stomatal conductance to various leaf water potential
    drop and thus, how likely the plant is to cavitate for a specific
    value of P.

    """

    return np.maximum(cst.zero, np.exp(- (- P / b) ** c))


def k_regulate(kP0, fP0):

    """
    Maximum xylem hydraulic conductance can be defined as the logically
    equivalent to the conductance corresponding to the maximum water
    potential.

    Arguments:
    ----------
    K0: float
        hydraulic conductance [mmol.s-1.m-2.MPa-1] at the upper
        potential bound of the water potential range, i.e. f(Ps=0, b, c)

    fP0: float
        current upper bound of the water potential range [MPa]

    Returns:
    --------
    The maximum hydraulic conductance [mmol.s-1.m-2.MPa-1] of the plant,
    given the current upper bound of the water potential range.

    """

    return kP0 * fP0


def Weibull_params(p):

    """
    Calculates the two Weibull parameters b and c, using two different
    leaf water potentials that each cause a decline of x% in hydraulic
    conductance. The function looks for two different x% values of P in
    the parameters, with x1 < x2. Such that:
        at x1% ln(1) - ln(1 - x1) = -ln(1 - x1) = ln(2)
        and x2% ln(1) - ln(1 - x2) = -ln(1 - x2) = -ln(1 / 50) = ln(50)

    Arguments:
    ----------
    Px1: float
        leaf water potential [MPa] causing a x1% decline in hydraulic
        conductance, with x1 < x2

    Px2: float
        leaf water potential [MPa] causing a x2% decline in hydraulic
        conductance, with x2 > x1

    Returns:
    --------
    b: float
        one of two Weibull parameters [MPa], it is P at k / kmax = 0.37

    c: float
        one of two Weibull parameters [unitless], it controls whether
        the curve has a threshold sigmoidal form or non-threshold
        sigmoidal form

    """

    # empty array to save the two x%
    x = []

    # find the x% and associated Px values
    for params in list(p.index):

        if str(params)[0] == 'P':

            try:
                if (int(str(params)[1:3])) and (len(str(params)) == 3):
                    x += [float(str(params)[1:3])]

            except ValueError:
                pass

    x1 = min(x)  # x1 < x2, at this stage, still a %
    x2 = max(x)
    Px1 = p['P%d' % (int(x1))]  # get the % value in the pandas series
    Px2 = p['P%d' % (int(x2))]
    x1 /= 100.  # normalise between 0-1
    x2 /= 100.

    try:  # c is derived from both expressions of b
        c = np.log(np.log(1. - x1) / np.log(1. - x2)) / (np.log(Px1) -
                                                         np.log(Px2))

    except ValueError:
        c = np.log(np.log(1. - x2) / np.log(1. - x1)) / (np.log(Px2) -
                                                         np.log(Px1))

    # b = Px2 / ((- np.log(1 - x2)) ** (1. / c)), or
    b = Px1 / ((- np.log(1 - x1)) ** (1. / c))

    return b, c


def transpiration(P, kmax, b, c):

    """
    Calculates the transpiration at steady-state by integrating over all
    vulnerability curves in the soil-plant system, i.e. over the full
    range of water potentials P. The err associated with the quad
    integral function is relatively small and this method is equivalent
    to using a gamma function.

    Arguments:
    ----------
    P: array or sympy symbol
        leaf water potential [MPa], either an array of values from the
        soil water potential Ps to the critical water potential Pcrit
        for which cavitation of the xylem occurs, or a sympy symbol to
        use in a symbolic expression

    kmax: float
        maximum hydraulic conductance [mmol s-1 m-2 MPa-1] of the plant

    b: float
        one of two Weibull parameters [MPa], it is P at k / kmax = 0.37

    c: float
        one of two Weibull parameters [unitless], it controls whether
        the curve has a threshold sigmoidal form or non-threshold
        sigmoidal form

    Returns:
    --------
    trans: array or sympy expression
        transpiration rate [mol.m-2.s-1]

    """

    trans = np.empty_like(P)  # empty numpy array of right length

    for i in range(len(P)):  # at Ps, trans=0; at Pcrit, trans=transcrit

        trans[i], err = quad(f, P[i], P[0], args=(b, c))

    trans[trans > cst.zero] *= kmax * conv.FROM_MILI  # mol.s-1.m-2

    try:
        return np.maximum(cst.zero, trans)

    except TypeError:
        return trans


def hydraulics(p, res='low'):

    """
    Calculates the hydraulics used to solve the leaf energy balance.
    These are the pressure drop range from Ps to Pcrit and the
    transpiration. Pcrit is defined as the critical point related to the
    critical conductance below which embolism occurs:
    Pcrit = - b * np.log(kmax / kcrit) ** (1. / c)
        with kcrit a percentage of kmax: kcrit = p.ratiocrit * kmax
    hence:
        Pcrit = - b * np.log(kmax / (p.ratiocrit * kmax)) ** (1. / c)

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    res: string
        either 'low' (default), 'med', or 'high' for the transpiration
        stream

    Returns:
    --------
    P: array or sympy symbol
        leaf water potential [MPa], either an array of values from the
        soil water potential Ps to the critical water potential Pcrit
        for which cavitation of the xylem occurs, or a sympy symbol to
        use in a symbolic expression

    trans: array or sympy expression
        transpiration rate [mol m-2 s-1]

    """

    # two Weibull parameters setting the shape of the vuln curves
    b, c = Weibull_params(p)  # MPa, unitless

    # get the leaf water potentials P
    Pcrit = - b * np.log(1. / p.ratiocrit) ** (1. / c)  # MPa

    if p.Ps <= Pcrit:  # plants die, the optimization cannot happen
        raise IndexError('plants are dead, no optimization')

    if res == 'low':
        P = np.linspace(p.Ps, Pcrit, 200)  # MPa, between Ps and Pcrit

    if res == 'med':
        P = np.linspace(p.Ps, Pcrit, 400)  # MPa

    if res == 'high':
        P = np.linspace(p.Ps, Pcrit, 600)  # MPa

    # integrate the vulnerability functions on P to get E
    trans = transpiration(P, p.kmax, b, c)  # mol m-2 s-1

    return P, trans
