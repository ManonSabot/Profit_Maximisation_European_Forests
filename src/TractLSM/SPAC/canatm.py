# -*- coding: utf-8 -*-

"""
Functions related to atmospheric processes relevant at the canopy level.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Abramowitz, G., Pouyanné, L., & Ajami, H. (2012). On the information
  content of surface meteorology for downward atmospheric long‐wave
  radiation synthesis. Geophysical Research Letters, 39(4).
* Dai et al. (2004) Journal of Climate, 17, 2281-2299.
* De Pury & Farquhar (1997). PCE, 20, 537-557.
* Jarvis, P. G., & McNaughton, K. G. (1986). Stomatal control of
  transpiration: scaling up from leaf to region. In Advances in
  ecological research (Vol. 15, pp. 1-49). Academic Press.
* Lloyd et al. (2010). Biogeosciences, 7, 1833–1859
* Norman, J. M., & Campbell, G. S. (1998). An introduction to
  environmental biophysics. Springer, New York.
* Medlyn et al. (2007). Linking leaf and tree water use with an
  individual-tree model. Tree Physiology, 27(12), 1687-1699.
* Monteith, J. L., & Unsworth, M. H. (1990). Principles of environmental
  physics. Arnold. SE, London, UK.
* Spitters et al. (1986). Separating the diffuse and direct component of
  global radiation and its implications for modeling canopy
  photosynthesis. Part I. Components of incoming radiation. Agricultural
  Forest Meteorol., 38:217-229.
* Wang and Leuning (1998) AFm, 91, 89-111.

"""

__title__ = "Canopy atmospheric processes"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (10.07.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import numpy as np  # array manipulations, math operators
from sympy import exp  # express exponential functions symbolically

# own modules
from TractLSM import conv, cst  # unit converter & general constants


# ======================================================================

def vpsat(T):

    """
    Calculates the saturation vapour pressure at a specific temperature
    T as given in Monteith & Unsworth, 1990.

    Arguments:
    ----------
    T: array or float
        temperature [degC]

    Returns:
    --------
    The saturation vapour pressure [kPa] at T.

    """

    try:
        return 0.61078 * np.exp(17.27 * T / (T + 237.3))

    except AttributeError:
        return 0.61078 * exp(17.27 * T / (T + 237.3))


def slope_vpsat(p):

    """
    Calculates the slope of saturation vapour pressure of water at air
    temperature.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    The slope of saturation vapour pressure of water [kPa degC-1].

    """

    return (vpsat(p.Tair + 0.1) - vpsat(p.Tair)) / 0.1


def LH_water_vapour(p):

    """
    Calculates the latent heat of water vapor at air temperature as
    given by eq A5 of Medlyn et al., 2007.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    The latent heat of water vapor [J mol-1].

    """

    return (cst.LH2O - 2.365e3 * p.Tair) * cst.MH2O * conv.FROM_MILI


def psychometric(p):

    """
    Calculates the atmospheric psychrometric constant.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    The psychrometric constant [kPa degC-1].

    """

    return cst.Cp * p.Patm / LH_water_vapour(p)


def emissivity(p):

    """
    Calculates the emissivity of the atmosphere by deriving it from the
    empirical long-wave down estimate proposed by Abramowitz et al.
    (2012): LWdown = 0.031 * ea + 2.84 * T - 522.5 (W m-2)

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    The apparent emissivity at air temperature [unitless].

    """

    # unit conversion
    TairK = p.Tair + conv.C_2_K  # degK

    # actual vapour pressure
    ea = (vpsat(p.Tair) - p.VPD) * conv.MILI  # Pa

    return (0.031 * ea + 2.84 * TairK - 522.5) / (cst.sigma * TairK ** 4.)


def net_radiation(p, surface='veg'):

    """
    Calculates net isothermal radiation, i.e. the net radiation that
    would be recieved if object and air temperature were the same
    (Eq 11.14 of Campbell and Norman, 1998.). Having a dependency on
    VPD through the emissivity partly accounts for clouds.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    surface: string
        'veg' or 'soil'

    Returns:
    --------
    The net isothermal radiation [W m-2].

    """

    # unit conversions
    TairK = p.Tair + conv.C_2_K  # degK

    # type of surface
    if surface == 'soil':
        emis_surf = p.eps_s
        albedo = p.albedo_s

    else:
        emis_surf = p.eps_l
        albedo = p.albedo_l

    # incoming short and long wave radiation
    Rsw = (1. - albedo) * p.PPFD * conv.PAR_2_SW  # W m-2
    Rlw = emissivity(p) * cst.sigma * TairK ** 4.  # W m-2
    Rabs = Rsw + Rlw  # absorbed radiation, W m-2

    return Rabs - emis_surf * cst.sigma * TairK ** 4.


def estimate_clearness(p):

    """
    Estimate atmospheric transmisivity - the amount of diffuse radiation
    is a function of the amount of haze and/or clouds in the sky.
    Estimate a proxy for this, i.e. the ratio between global solar
    radiation on a horizontal surface at the ground and the
    extraterrestrial solar radiation.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    Atmospheric transmisivity (0.-1.)

    """

    # rad incident outside earth atmosphere, eq 1 Spitters et al. (1986)
    Sxt = cst.S0 * p.coszen * (1. + 0.033 * np.cos(2. * np.pi * p.doy /
                                                   365.25))

    # catch possible divide by zero when zenith = 90.
    if Sxt <= 0.0:
        tau = 0.0

    if (p.coszen > 1.e-10) and (p.PPFD * conv.PAR_2_SW > 10.):
        tau = p.PPFD * conv.PAR_2_SW / Sxt

    else:
        tau = 0.

    return np.maximum(0., np.minimum(1., tau))


def spitters(p):

    """
    Algorithm to estimate the diffuse component from the measured
    irradiance as given in eq 20a-d of Spitters et al., 1986.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    direct_frac: float
        the direct component of incoming radiation

    diffuse_frac: float
        the diffuse component of incoming radiation

    """

    # atmospheric transmisivity
    tau = estimate_clearness(p)

    # Spitters formula
    R = 0.847 + p.coszen * (1.04 * p.coszen - 1.61)

    if p.coszen <= 1.e-2:
        diffuse_frac = 1.

    elif (tau > 0.22) and (tau <= 0.35):
        diffuse_frac = np.maximum(0., np.minimum(1., 1. -
                                                 6.4 * (tau - 0.22) ** 2.))

    elif (tau > 0.35) and (tau <= (1.47 - R) / 1.66):
        diffuse_frac = np.maximum(0., np.minimum(1., 1.4728 - 1.66 * tau))

    else:
        diffuse_frac = np.maximum(0., np.minimum(1., R))

    diffuse_frac = np.maximum(0., np.minimum(1., diffuse_frac))

    return 1. - diffuse_frac, diffuse_frac


def psi_func(lai, z):

    """
    B5 function from Wang and Leuning (1998) which integrates a property
    over the canopy space. To avoid underflow, z * lai cannot be greater
    than 30!

    Arguments:
    ----------
    lai: float
        leaf area index [m2 m-2]

    z: float
        canopy property

    Returns:
    --------
    LAI upscaling coefficient

    """

    return (1. - np.exp(-np.minimum(z * lai, 30.))) / z


def scaling_2_leaves(p, kb, kd):

    """
    Calculate scaling factors from a big leaf canopy to a two leaves
    canopy, following Leuning et al. (1995), De Pury & Farquhar (1997),
    Wang and Leuning (1998).

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    kb: float
        black leaves direct extinction coef.

    kd: float
        black leaves diffuse extinction coef.

    Returns:
    --------
    fLAI: array
        LAI of the sunlit and shaded leaves [m2 m-2]

    fscale2can: array
        scaling term to go from big leaf to canopy [m2 m-2]

    fgradis: array
        scaling term for leaf radiative conductance [mol m-2 s-1]

    """

    # unit conversions
    TairK = p.Tair + conv.C_2_K  # degK
    cmolar = p.Patm * conv.MILI / (cst.R * TairK)  # air molar density

    # declare empty arrays
    fLAI = np.zeros(2)  # sunlit, shaded
    fscale2can = np.zeros(2)
    fgradis = np.zeros(2)

    if p.LAI > cst.LAI_thresh:  # vegetated
        transb = np.exp(-(np.minimum(kb * p.LAI, 30.)))  # no underflow
        transd = np.exp(-kd * p.LAI)

        # radiative conductance following eq D7 in Leuning et al (1995)
        a1 = 4. * p.eps_l / cst.Cp
        a2 = cst.sigma * kd * TairK ** 3
        a3 = (1. - transb * transd) / (kb + kd) + (transd - transb) / (kb - kd)
        fgradis[0] = np.maximum(1.e-3, a1 * a2 * a3)

        a1 = 2. * a1
        a3 = (1. - transd) / kd - fgradis[0]
        fgradis[1] = np.maximum(1.e-3, a1 * a2 * a3)

        # from m s-1 to mol m-2 s-1
        fgradis *= cmolar

    # sunlit & shaded LAI fracs of can (eq 18; De Pury & Farquhar, 1997)
    if (p.LAI > cst.LAI_thresh) and (p.PPFD * conv.PAR_2_SW > cst.RAD_thresh):
        fLAI[0] = (1. - np.exp(-kb * p.LAI)) / kb

    else:
        fLAI[0] = 0.

    fLAI[1] = p.LAI - fLAI[0]

    # from big sunlit & shaded leaves to can (Wang and Leuning, 1998)
    fscale2can[0] = (1. - np.exp(-(kb + p.kn) * p.LAI)) / (kb + p.kn)
    fscale2can[1] = (1. - np.exp(-p.kn * p.LAI)) / p.kn - fscale2can[0]

    return fLAI, fscale2can, fgradis


def absorbed_radiation_2_leaves(p):

    """
    Calculate absorded irradiance of sunlit and shaded fractions of the
    canopy. The Ross-Goudriaan algorithm is used to get the ratio of the
    projected area of leaves in the direction perpendicular to the
    direction of incident solar radiation and the actual leaf area
    (Sellers (1985) eq. 13; CABLE code, Kowalczyk (2006) eq. 28/29).

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    fRcan: array
        Rnet absorbed by each leaf in the canopy [W m-2]

    fPPFD: array
        PPFD of the sunlit and shaded leaves [umol m-2 s-1]

    fLAI: array
        LAI of the sunlit and shaded leaves [m2 m-2]

    fscale2can: array
        scaling term to go from big leaf to canopy [m2 m-2]

    fgradis: array
        scaling term for leaf radiative conductance [mol m-2 s-1]

    """

    # unit conversions
    TairK = p.Tair + conv.C_2_K  # degK

    # declare empty arrays
    albedo_soil = np.zeros(2)  # sunlit, shaded
    spec = np.zeros(3)  # visible, NIR, LW
    fRcan = np.zeros((2, 3))  # both leaves and all light ranges

    # set values for use in the function
    cos3 = np.array([np.cos(np.deg2rad(15.)), np.cos(np.deg2rad(45.)),
                     np.cos(np.deg2rad(75.))])  # 15, 45, 75 deg angles
    gauss_w = np.array([0.308, 0.514, 0.178])  # Gaussian integ. weights
    tau_l = np.array([p.tau_l, 0.3])   # leaf transmis. VIS, NIR
    refl_l = np.array([p.albedo_l, 0.3])  # leaf reflect. VIS, NIR

    # direct/beam and diffuse fracs
    direct_frac, diffuse_frac = spitters(p)

    # leaf angle parameters
    xphi1 = 0.5 - p.chi_l * (0.633 + 0.33 * p.chi_l)
    xphi2 = 0.877 * (1. - 2. * xphi1)

    # Ross-Goudriaan algo (approximated as eq 28; Kowalcyk et al., 2006)
    RG = xphi1 + xphi2 * p.coszen

    # direct beam extinction coef. (eq 26; Kowalcyk et al., 2006)
    if (p.LAI > cst.LAI_thresh) and (direct_frac > cst.RAD_thresh):
        kb = RG / p.coszen

    else:  # i.e. bare soil
        kb = 0.5

    # diffuse extinction coef. (eq 27; Kowalcyk et al., 2006)
    if p.LAI > cst.LAI_thresh:
        kbx = (xphi1 + xphi2 * cos3) / cos3  # approx integration of kb
        kd = -np.log(np.nansum(gauss_w * np.exp(-kbx * p.LAI))) / p.LAI

    else:  # i.e. bare soil
        kd = 0.7

    if np.abs(kb - kd) < cst.RAD_thresh:
        kb = kd + cst.RAD_thresh

    if direct_frac < cst.RAD_thresh:
        kb = 1.e5

    # scaling coefficients, from big leaf canopy to two leaf canopy
    fLAI, fscale2can, fgradis = scaling_2_leaves(p, kb, kd)

    # reflect black horiz leaves (eq 6.19; Goudriaan & van Laar, 1994)
    spec[:2] = (1. - tau_l - refl_l) ** 0.5
    spec[2] = 1.
    crefb = (1. - spec) / (1. + spec)

    # reflect of diffuse radiation for black leaves
    if p.LAI > cst.LAI_thresh:
        crefd = 2. * crefb * (gauss_w[0] * kbx[0] / (kbx[0] + kd) +
                              gauss_w[1] * kbx[1] / (kbx[1] + kd) +
                              gauss_w[2] * kbx[2] / (kbx[2] + kd))

    else:  # no veg, no reflection
        crefd = np.zeros(3)

    # calculate soil albedo under leaves (ignoring snow)
    if p.albedo_s <= 0.14:
        fact = 0.5

    elif (p.albedo_s > 0.14) and (p.albedo_s <= 0.2):
        fact = 0.62

    else:
        fact = 0.68

    albedo_soil[1] = 2. * p.albedo_s / (1. + fact)  # shaded
    albedo_soil[0] = fact * albedo_soil[1]  # sunlit

    # update extinct coefs and fract transmittance (not black leaves)
    kd_m = kd * spec[:2]  # modified kd for leaf scattering (eq 6.20)
    ctransd = np.exp(-kd_m * p.LAI)  # diffuse transmittance fraction

    # calculate effective canopy-soil diffuse reflectance fraction
    if p.LAI > cst.LAI_thresh:
        ecrefd = crefd[:2] + (albedo_soil - crefd[:2]) * ctransd ** 2.

    else:  # i.e. bare soil
        ecrefd = albedo_soil

    # modified kb for leaf scattering (eq 6.20)
    if (p.LAI > cst.LAI_thresh) and (p.PPFD * conv.PAR_2_SW > cst.RAD_thresh):
        kb_m = kb * spec[:2]

    else:
        kb_m = 1.e-9

    # canopy direct transmittance fraction
    ctransb = np.exp(-np.minimum(kb_m * p.LAI, 30.))  # no underflow

    # modified canopy beam reflectance (6.21)
    crefb_m = 2. * kb / (kb + kd) * crefb[:2]

    # calculate effective canopy-soil direct reflectance fraction
    ecrefb = crefb_m + (albedo_soil - crefb_m) * ctransb ** 2.

    if (p.LAI > cst.LAI_thresh) and (p.PPFD * conv.PAR_2_SW > cst.RAD_thresh):
        Id = p.PPFD * conv.UMOL_2_J * diffuse_frac  # W m-2
        Ib = p.PPFD * conv.UMOL_2_J * direct_frac  # W m-2

        # rad absorbed by sunlit leaf(eq B3b; Wang and Leuning, 1998)
        a1 = Id * (1. - ecrefd) * kd_m
        a2 = psi_func(p.LAI, kd_m + kb)
        a3 = Ib * (1. - ecrefb) * kb_m
        a4 = psi_func(p.LAI, kb_m + kb)
        a5 = Ib * (1. - tau_l - refl_l) * kb
        a6 = psi_func(p.LAI, kb) - psi_func(p.LAI, 2. * kb)
        fRcan[0, :2] = a1 * a2 + a3 * a4 + a5 * a6

        # rad absorbed by shaded leaf, (eq B4; Wang and Leuning, 1998)
        a2 = psi_func(p.LAI, kd_m) - a2
        a4 = psi_func(p.LAI, kb_m) - a4
        fRcan[1, :2] = a1 * a2 + a3 * a4 - a5 * a6

    # longwave rad absorbed by sunlit, isothermal conditions (eq B18)
    a1 = -kd * cst.sigma * TairK ** 4.
    a2 = p.eps_l * (1. - emissivity(p))
    a3 = psi_func(p.LAI, kb + kd)
    a4 = 1. - p.eps_s
    a5 = (p.eps_l - emissivity(p))
    a6 = psi_func(p.LAI, 2. * kd) * psi_func(p.LAI, kb - kd)
    fRcan[0, 2] = a1 * (a2 * a3 + a4 * a5 * a6)

    # longwave rad absorbed by shaded, isothermal conditions (eq B19)
    a3 = psi_func(p.LAI, kd)
    a6 = np.exp(-kd * p.LAI) * a3
    fRcan[1, 2] = a1 * (a2 * a3 - a4 * a5 * a6) - fRcan[0, 2]

    # PPFD of the sunlit and shaded leaves in umol m-2 s-1
    fPPFD = fRcan[:, 0] * conv.J_2_UMOL  # umol m-2 s-1

    # total energy absorbed by canopy, summing VIS, NIR and LW
    fRcan = np.nansum(fRcan, axis=1)  # net isothermal radiation of can

    return fRcan, fPPFD, fLAI, fscale2can, fgradis
