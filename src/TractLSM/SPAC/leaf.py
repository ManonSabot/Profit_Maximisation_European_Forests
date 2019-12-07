# -*- coding: utf-8 -*-

"""
Functions related to leaf canopy processes: used to calculate
conductances, leaf temperature, and photosynthesis.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Collatz et al. (1991). Regulation of stomatal conductance and
  transpiration: a physiological model of canopy processes. Agric. For.
  Meteorol, 54, 107-136.
* De Pury, D. G. G., & Farquhar, G. D. (1997). Simple scaling of
  photosynthesis from leaves to canopies without the errors of big‐leaf
  models. Plant, Cell & Environment, 20(5), 537-557.
* Farquhar, G. D., von Caemmerer, S. V., & Berry, J. A. (1980). A
  biochemical model of photosynthetic CO2 assimilation in leaves of C3
  species. Planta, 149(1), 78-90.
* Jones, H. G. (2013). Plants and microclimate: a quantitative approach
  to environmental plant physiology. Cambridge university press.
* Kirschbaum, M. U. F., & Farquhar, G. D. (1984). Temperature dependence
  of whole-leaf photosynthesis in Eucalyptus pauciflora Sieb. ex Spreng.
  Functional Plant Biology, 11(6), 519-538.
* Leuning, R. (1990). Modelling stomatal behaviour and photosynthesis of
  Eucalyptus grandis. Functional Plant Biology, 17(2), 159-175.
* Norman, J. M., & Campbell, G. S. (1998). An introduction to
  environmental biophysics. Springer, New York.
* Medlyn et al. (2002). Temperature response of parameters of a
  biochemically based model of photosynthesis. II. A review of
  experimental data. Plant, Cell & Environment, 25(9), 1167-1179.
* Medlyn et al. (2007). Linking leaf and tree water use with an
  individual-tree model. Tree Physiology, 27(12), 1687-1699.
* Monteith, J. L., & Unsworth, M. H. (1990). Principles of environmental
  physics. Arnold. SE, London, UK.
* Slatyer, R. O. (1967). Plant-water relationships. Academic Press; New
  York; San Frncisco; London.
* Tjoelker et al. (2001). Modelling respiration of vegetation: evidence
  for a general temperature‐dependent Q10. Global Change Biology, 7(2),
  223-230.

"""

__title__ = "Leaf level photosynthetic processes"
__author__ = "Manon E. B. Sabot"
__version__ = "2.0 (10.07.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import numpy as np  # array manipulations, math operators
from sympy import exp, Min  # express functions symbolically

# own modules
from TractLSM import conv, cst  # unit converter & general constants
from TractLSM.SPAC import vpsat, slope_vpsat, LH_water_vapour


# ======================================================================

def conductances(p, Tleaf=None, gs=None):

    """
    Both forced and free convection (gHa) contribute to exchange of heat
    and mass through leaf boundary layers at the wind speeds typically
    encountered within plant canopies (< 0-5 m.s-1). Total leaf
    conductance to heat and total leaf conductance to water vapour (or
    simply boundary layer conductance to water vapour) are needed for
    the energy balance. The leaf LAI is used in gHf to adjust for the
    2-leaf model.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Tleaf: array or float or sympy expression
        leaf temperature [degC]

    gs: float
        stomatal conductance to water vapour [mol m-2 s-1]

    Returns:
    --------
    gw: float
            total leaf conductance to water vapour [mol m-2 s-1]

    gH: float
        total leaf conductance to heat [mol m-2 s-1]

    gb: float
        boundary layer conductance to water vapour [mol m-2 s-1]

    gr: float
        radiative conductance [mol m-2 s-1]

    """

    # unit conversions
    TairK = p.Tair + conv.C_2_K  # degK
    cmolar = p.Patm * conv.MILI / (cst.R * TairK)  # air molar density

    # Sutherland Eq for dynamic viscosity
    mu = 1.458e-6 * TairK ** 1.5 / (TairK + 110.4)  # Pa s

    # kinematic viscosity
    nu = mu * cst.R * TairK / (p.Patm * cst.Mair)  # m2 s-1
    prandtl = nu / cst.DH  # unitless

    # boundary layer cond to forced convect. (Campbell & Norman, 1998)
    d = 0.72 * p.max_leaf_width  # leaf width, m
    reynolds = p.u * d / nu  # unitless
    gHa = (0.664 * cmolar * cst.DH * reynolds ** 0.5 * prandtl ** (1. / 3.)
           / d)

    if Tleaf is not None:
        grashof = ((cst.g0 * (1. / TairK) * abs(Tleaf - p.Tair) * d ** 3.) /
                   nu ** 2.)  # unitless

        # boundary layer cond to free convect. (Campbell & Norman, 1998)
        gHf = (p.LAI * (0.54 * cmolar * cst.DH * (grashof * prandtl) ** 0.25) /
               d)

    else:
        gHf = 1.e-3  # default minimum value, from CABLE

    try:
        gHa = np.maximum(1.e-3, gHa + gHf)  # mol m-2 s-1

    except TypeError:  # if used as expression, np.maximum won't work
        gHa = gHa + gHf  # mol m-2 s-1

    # boundary layer conductance to water vapour
    if not np.isclose(gHa, cst.zero, rtol=cst.zero, atol=cst.zero):
        gb = np.maximum(cst.zero, gHa * conv.GbvGbh)  # mol m-2 s-1

    else:
        gb = cst.zero  # mol m-2 s-1

    # radiative conductance (Jones, 2013; Medlyn et al., 2007)
    gr = np.maximum(cst.zero, 4. * p.eps_l * cst.sigma * TairK ** 3. / cst.Cp)

    # total two-sided leaf conductance to heat (Medlyn et al., 2007)
    if (np.isclose(gHa, cst.zero, rtol=cst.zero, atol=cst.zero) and
       np.isclose(gr, cst.zero, rtol=cst.zero, atol=cst.zero)):
        gH = cst.zero

    else:
        gH = 2. * (gHa + gr)  # hypostomatous leaf (2-sided)

    if gs is None:

        return gH, gb, gr

    if gs is not None:

        # total leaf conductance to water vapour (Medlyn et al., 2007)
        if (gs > -cst.zero) and not np.isclose(gb, cst.zero, rtol=cst.zero,
           atol=cst.zero) and not np.isclose(gs, cst.zero, rtol=cst.zero,
           atol=cst.zero):
            gw = np.maximum(cst.zero, (gb * gs) / (gb + gs))

        else:
            gw = cst.zero

        return gw, gH, gb, gr


def leaf_temperature(p, trans, Tleaf=None, gradis=False):

    """
    Calculates the leaf temperature for each supply function, i.e. over
    the transpiration stream. A factor 2 is introduced in the
    denominator because gHa and gr because leaves are two-sided.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    trans: array or float or sympy symbol
        transpiration rate [mol m-2 s-1]

    Tleaf: array or float or sympy expression
        leaf temperature [degC]

    gradis: array or float
        scaling term for leaf radiative conductance [mol m-2 s-1]

    Returns:
    --------
    Tleaf: array or float or sympy expression
        leaf temperature [degC]

    gb: float
        boundary layer conductance to water vapour [mol m-2 s-1]

    """

    # unit conversion
    TairK = p.Tair + conv.C_2_K  # degK

    # get conductances
    gH, gb, __ = conductances(p, Tleaf=Tleaf)  # mol m-2 s-1

    # latent heat of water vapor
    Lambda = LH_water_vapour(p)  # J mol-1

    # slope of saturation vapour pressure of water
    slp = slope_vpsat(p)  # kPa degK-1

    # canopy/leaf sensible heat flux
    H = p.Rnet - Lambda * trans  # W m-2

    if gradis:  # update canopy/leaf sensible heat flux
        TleafK = Tleaf + conv.C_2_K  # degK
        H -= (TleafK - TairK) * cst.Cp * gH * p.gradis / (gH + p.gradis)

    # simplified Tleaf (gb for gw), eq 14.6 of Campbell & Norman, 1998
    if np.isclose(abs(p.Tair), 0., rtol=cst.zero, atol=cst.zero):
        Tleaf = (p.Tair + H / (cst.Cp * gH * TairK / cst.zero +
                               Lambda * slp * gb / p.Patm))  # degC
    else:
        Tleaf = (p.Tair + H / (cst.Cp * gH * TairK / p.Tair +
                               Lambda * slp * gb / p.Patm))  # degC

    if gradis:  # update net canopy/leaf radiation
        TleafK = Tleaf + conv.C_2_K  # degK
        p.Rnet -= (TleafK - TairK) * cst.Cp * p.gradis  # W m-2

    return Tleaf, gb


def leaf_energy_balance(p, trans, Tleaf=None):

    """
    Calculates the CO2 diffusive conductance of leaf gc using the
    saturation vapour pressure deficit of water (vpsat) and the stomatal
    conductance to water vapour (gs) using the differential between the
    water content of saturated air at leaf temperature and at air
    temperature.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    trans: array or float or sympy symbol
        transpiration rate [mol m-2 s-1]

    Tleaf: array or float or sympy expression
        leaf temperature [degC]

    Returns:
    --------
    gc: array or sympy expression
        leaf CO2 diffusive leaf conductance [mol s-1 m-2]

    gs: array
        stomatal conductance to water vapour [mol m-2 s-1]

    gb: float
        boundary layer conductance to water vapour [mol m-2 s-1]

    ww: array
        plant-air saturated H2O content differential
        [moles(H2O) mole-1(air)]

    """

    if Tleaf is None:
        Tleaf, gb = leaf_temperature(p, trans)

    else:
        __, gb = leaf_temperature(p, trans, Tleaf=Tleaf)

    # saturation vapour pressure of water at T
    esat_l = vpsat(Tleaf)  # kPa
    esat_a = vpsat(p.Tair)  # kPa

    # get the leaf CO2 diff. conductance
    Dleaf = (esat_l - (esat_a - p.VPD))  # leaf-air vpd, kPa
    gw = p.Patm * trans / Dleaf  # leaf vapour diff. cond, mol s-1 m-2

    try:
        gc = np.maximum(cst.zero, gw * conv.GcvGw)  # mol s-1 m-2
        gw[np.isclose(trans, cst.zero, rtol=cst.zero, atol=cst.zero)] = \
            cst.zero
        gc[np.isclose(gw, cst.zero, rtol=cst.zero, atol=cst.zero)] = cst.zero

    except TypeError:  # if used as expression, np.maximum won't work
        gc = gw * conv.GcvGw  # mol s-1 m-2

    # partial pressure of dry air
    pa = conv.FROM_MILI * p.Patm * cst.Mair * cst.Rd / cst.R  # kPa

    # partial pressure ot water vap
    pw = p.Patm - pa  # kPa

    # plant-air H2O vap diff (Slatyer, 1967), moles(H2O) mole-1(air)
    ww = (cst.MH2O / cst.Mair) * (esat_l / (pa - esat_l) - pw / (pa - pw))

    # gs, stomatal conductance to water vapour
    gs = gb * gw / (gb - gw)  # mol s-1 m-2

    try:
        gs[gs < 0.] = cst.zero
        gs[np.isclose(trans, cst.zero, rtol=cst.zero, atol=cst.zero)] = \
            cst.zero

    except TypeError:  # if used as expression
        pass

    return gc, gs, gb, ww


def arrhen(v25, Ea, Tref, Tleaf, deltaS=None, Hd=None):

    """
    Calculates the temperature dependence of a kinetic variable using an
    Arrhenius function which transforms the variable at 25 degC given
    its energy of activation and the leaf temperature (Medlyn et al.,
    2002). Providing deltaS and Hd returns a peaked Arrhenius function
    which accounts for the rate of inhibition at higher temperatures.

    Arguments:
    ----------
    v25: float
        kinetic variable at 25 degC [varies]

    Ea: float
        energy of activation of the variable [J mol-1]

    Tleaf: array or float or sympy expression
        leaf temperature [degC]

    deltaS: float
        entropy factor [J mol-1 K-1]

    Hd: float
        rate of decrease about the optimum temperature [J mol-1]

    Returns:
    --------
    The temperature-dependent kinetic variable [varies].

    """

    # unit conversion
    Tl_K = Tleaf + conv.C_2_K  # degK

    try:
        arrhenius = v25 * np.exp(Ea * (Tl_K - Tref) / (Tref * cst.R * Tl_K))

    except AttributeError:
        arrhenius = v25 * exp(Ea * (Tl_K - Tref) / (Tref * cst.R * Tl_K))

    if (deltaS is None) or (Hd is None):
        return arrhenius

    else:
        try:
            arg2 = 1. + np.exp((deltaS * Tref - Hd) / (cst.R * Tref))
            arg3 = 1. + np.exp((deltaS * Tl_K - Hd) / (cst.R * Tl_K))

        except AttributeError:
            arg2 = 1. + exp((deltaS * Tref - Hd) / (cst.R * Tref))
            arg3 = 1. + exp((deltaS * Tl_K - Hd) / (cst.R * Tl_K))

        return arrhenius * arg2 / arg3


def adjust_low_T(var, Tleaf, lower_bound=0., upper_bound=10.):

    """
    Function linearly forcing a variable to zero at low temperature

    Arguments:
    ----------
    var: float or array
        kinetic variable [varies]

    Tleaf: array or float or sympy expression
        leaf temperature [degC]

    lower_bound: float
        lowest possible leaf temperature [degC]

    upper_bound: float
        upper "lower" leaf temperature [degC]

    Returns:
    --------
    The temperature-dependent kinetic variable [varies].

    """

    if 'float' in str(type(Tleaf)):
        if Tleaf < lower_bound:
            var = 0.

        elif Tleaf < upper_bound:
            var *= (Tleaf - lower_bound) / (upper_bound - lower_bound)

    else:
        if np.any(Tleaf < lower_bound):
            var[Tleaf < lower_bound] = 0.

        if np.any(Tleaf < upper_bound):
            low = np.logical_and(Tleaf < upper_bound, Tleaf > lower_bound)

            try:
                var[np.where(low)] *= ((Tleaf[np.where(low)] - lower_bound) /
                                       (upper_bound - lower_bound))

            except TypeError:
                var[low] *= ((Tleaf[low] - lower_bound) /
                             (upper_bound - lower_bound))

    return var


def foliar_resp(p, Tleaf):

    """
    Calculates the foliar respiration of plants given a reference
    temperature and a reference respiration rate, as proposed by
    Tjoelker et al., 2001.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Tleaf: array or float or sympy expression
        leaf temperature [degC]

    Returns:
    --------
    The leaf respiration rate [μmol m-2 s-1] at air temperature.

    """

    Q10 = 3.22 - 0.046 * Tleaf  # ratio of respiration T:(T - 10.)

    return p.Rlref * Q10 ** ((Tleaf - p.TRlref) / 10.)


def quad(a, b, c, large_root=True):

    """
    Calculates the square root given by the quadratic formula,
        with a, b, and c from ax2 + bx + c = 0.

    Arguments:
    ----------
    a, b, c: float or sympy symbolic expression
        coefficients of the equation to solve

    large_root: boolean
        if True, the largest root is returned

    Returns:
    --------
    Either one of the large or small roots given by the quadratic
    formula.

    """

    if large_root:
        return 0.5 * (-b + (b ** 2. - 4. * a * c) ** 0.5) / a

    else:
        return 0.5 * (-b - (b ** 2. - 4. * a * c) ** 0.5) / a


def quad_solve_Ci(Cs, gs_over_A, Rleaf, gamstar, v1, v2):

    """
    Solves for Ci starting from Cs, according to the standard quadratic
    way of solving for Ci as described in Leuning, 1990.

    Arguments:
    ----------
    Cs: float
        leaf surface CO2 concentration [Pa]

    gs_over_A: float
        gs/A as predicted by the USO (Medlyn, 2011) model

    Rleaf: float
        leaf dark respiration [μmol m-2 s-1]

    gamstar: float
        CO2 compensation point [Pa]

    v1: float
        Vmax or J

    v2: float
        Km or 2 * gamstar

    Returns:
    --------
    The intercellular CO2 concentration [Pa].

    """

    # unit conversions, from Pa to μmol mol-1
    Csi = Cs * conv.MILI * conv.FROM_kPa
    gammastar = gamstar * conv.MILI * conv.FROM_kPa
    V2 = v2 * conv.MILI * conv.FROM_kPa

    g0 = 1.e-9  # removing g0 introduces a solving error

    a = g0 + gs_over_A * (v1 - Rleaf)
    b = ((1. - Csi * gs_over_A) * (v1 - Rleaf) + g0 * (V2 - Csi) - gs_over_A *
         (v1 * gammastar + V2 * Rleaf))
    c = - ((1. - Csi * gs_over_A) * (v1 * gammastar + V2 * Rleaf) + g0 * V2 *
           Csi)

    ref_root = quad(a, b, c) * conv.ref_kPa * conv.FROM_MILI

    if (ref_root > Cs) or (ref_root < cst.zero):
        return quad(a, b, c, large_root=False) * conv.ref_kPa * conv.FROM_MILI

    else:
        return ref_root


def calc_photosynthesis(p, trans, Ci_s, photo, smooth=True, Tleaf=None,
                        Rleaf=None, gs_over_A=None):

    """
    Calculates the assimilation rate given the internal leaf CO2
    concentration following either the classic Farquhar photosynthesis
    model (with smoothed solve) or the Collatz model. The non-smoothed
    alternative would be An = min(Aj, Ac) - Rleaf.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    trans: array or float or sympy symbol
        transpiration rate [mol m-2 s-1]

    Ci_s: array or float or sympy symbol
        intercellular CO2 concentration [Pa] (the leaf surface CO2
        concentration [Pa] can be parsed instead with gs_over_A as well
        to solve for Ci)

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    smooth: boolean
        for the Sperry model to accurately solve for Ci, the transition
        point between Aj and Ac must be smoothed. True is the default

    Tleaf: float
        leaf temperature [degC]

    Rleaf: float
        leaf dark respiration [μmol m-2 s-1]

    gs_over_A: float
        gs/A as predicted by USO (Medlyn, 2011) model. Used for the
        quadratic solving of Ci

    Returns:
    --------
    A: array or float or sympy symbolic expression
        gross C assimilation rate [μmol m-2 s-1]

    Aj: array or float or sympy symbolic expression
        electron transport-limited photosynthesis rate [μmol m-2 s-1]

    Ac: array or float or sympy symbolic expression
        rubisco-limited photosynthesis rate [μmol m-2 s-1]

    """

    # initialise Ci
    Ci = Ci_s

    if Tleaf is None:
        Tleaf, __ = leaf_temperature(p, trans)

    if Rleaf is None:  # leaf dark respiration
        Rleaf = foliar_resp(p, Tleaf)  # μmol m-2 s-1

    # gamstar, Vmax, Kc and Ko are known at Tref, get their T dependency
    Tref = p.Tref + conv.C_2_K  # degk, Tref set to 25 degC

    # CO2 compensation point
    gamstar = arrhen(p.gamstar25, p.Egamstar, Tref, Tleaf)

    # max carbohylation rate, μmol m-2 s-1
    Vmax = arrhen(p.Vmax25, p.Ev, Tref, Tleaf, deltaS=p.deltaSv, Hd=p.Hdv)

    # max electron transport rate, μmol m-2 s-1
    Jmax = arrhen(p.JV * p.Vmax25, p.Ej, Tref, Tleaf, deltaS=p.deltaSj,
                  Hd=p.Hdj)

    if 'sympy' not in str(type(Tleaf)):  # adjust for low temperatures
        Vmax = adjust_low_T(Vmax, Tleaf)
        Jmax = adjust_low_T(Jmax, Tleaf)

    # Michaelis-Menten constants
    Kc = arrhen(p.Kc25, p.Ec, Tref, Tleaf)  # cst for carboxylation, Pa
    Ko = arrhen(p.Ko25, p.Eo, Tref, Tleaf)  # cst for oxygenation, kPa

    # at T = 0., the Ko can fail
    try:
        Ko = np.maximum(cst.zero, Ko)  # we don't want zeros in Km div

    except TypeError:
        pass

    # Michaelis-Menten constant for O2/CO2
    Km = Kc * (1. + p.O2 / Ko)

    # scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C
    Rleaf *= p.scale2can
    Vmax *= p.scale2can
    Jmax *= p.scale2can

    # rubisco-limited photosynthesis rate (De Pury & Farquhar, 1997)
    if gs_over_A is not None:
        Ci = quad_solve_Ci(Ci_s, gs_over_A, Rleaf, gamstar, Vmax, Km)
        Ci_c = Ci  # track Ci to know which one was used

    try:
        if (Ci <= cst.zero) or (Ci > Ci_s) or np.isnan(Ci):
            Ac = 0.
            Ci_c = 0.

        else:
            Ac = Vmax * (Ci - gamstar) / (Ci + Km)  # μmol m-2 s-1

    except (TypeError, ValueError):
        Ac = Vmax * (Ci - gamstar) / (Ci + Km)  # μmol m-2 s-1

    # electron transport-limited photosynthesis rate
    if photo == 'Farquhar':  # De Pury & Farquhar, 1997
        J = quad(p.c1, -((1. - p.tau_l - p.albedo_l) * p.alpha * p.PPFD +
                 Jmax), (1. - p.tau_l - p.albedo_l) * p.alpha * p.PPFD * Jmax,
                 large_root=False)

    else:  # Collatz et al. 1991
        J = (1. - p.tau_l - p.albedo_l) * p.alpha * p.PPFD

    # RuBP regeneration rate accounted for
    J *= 0.25  # μmol m-2 s-1

    if gs_over_A is not None:
        Ci = quad_solve_Ci(Ci_s, gs_over_A, Rleaf, gamstar, J, 2. * gamstar)

    Aj = J * (Ci - gamstar) / (Ci + 2. * gamstar)  # μmol m-2 s-1

    try:  # below light compensation point?
        if (Ci - gamstar <= cst.zero) or (Ci > Ci_s) or np.isnan(Ci):
            if gs_over_A is not None:
                Ci = Ci_s  # reinitialise Ci and solve using Cs
                Aj = J * (Ci - gamstar) / (Ci + 2. * gamstar)

            else:
                Aj = 0.

    except (TypeError, ValueError):
        pass

    if smooth:  # smooth transition point (Kirschbaum & Farquhar, 1984)
        if photo == 'Farquhar':
            An = ((Aj + Ac - ((Aj + Ac) ** 2. - 4. * p.c2 * Aj * Ac) ** 0.5) /
                  (2. * p.c2) - Rleaf)  # μmol m-2 s-1

        else:
            An = ((Aj + Ac - ((Aj + Ac) ** 2. - 4. * p.c4 * Aj * Ac) ** 0.5) /
                  (2. * p.c4) - Rleaf)  # μmol m-2 s-1

    else:  # non-smoothed transition pt between Aj and Ac
        if 'sympy' not in str(type(Tleaf)):
            An = min(Aj, Ac) - Rleaf

        else:
            An = Min(Aj, Ac) - Rleaf

    if gs_over_A is not None:
        Rublim = rubisco_limit(Aj, Ac)

        if Rublim == 'True':
            Ci = Ci_c  # the system is Rubisco-limited, this Ci was used

        return An, Aj, Ac, Ci

    else:
        return An, Aj, Ac


def rubisco_limit(Aj, Ac):

    """
    Tests whether the standard model for photosynthesis is rubisco
    limited or not, in which case it is limited by electron transport.

    Arguments:
    ----------
    Aj: float
        electron transport-limited photosynthesis rate [μmol m-2 s-1]

    Ac: float
        rubisco-limited photosynthesis rate [μmol m-2 s-1]

    Returns:
    --------
    'True' if the C assimilation is rubisco limited, 'False' otherwise.

    """

    if (np.minimum(Ac, Aj) > 0.) and np.isclose(np.minimum(Ac, Aj), Ac):
        return str(bool(1))

    elif (np.minimum(Ac, Aj) > 0.) and np.isclose(np.minimum(Ac, Aj), Aj):
        return str(bool(0))

    else:
        return 0.
