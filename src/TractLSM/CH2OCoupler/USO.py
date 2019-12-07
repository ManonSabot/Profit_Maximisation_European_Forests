# -*- coding: utf-8 -*-

"""
Imitate the way a LSM would solve for photosynthetic assimilation and
transpiration by following an iteration logic on the air temperature to
get the leaf temperature for which the Penman-Monteith energy balance
conditions are satisfied.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
----------
* Kowalczyk, E. A., Wang, Y. P., Law, R. M., Davies, H. L., McGregor,
  J. L., & Abramowitz, G. (2006). The CSIRO Atmosphere Biosphere Land
  Exchange (CABLE) model for use in climate models and as an offline
  model. CSIRO Marine and Atmospheric Research Paper, 13, 42.
* Medlyn, B. E., Duursma, R. A., Eamus, D., Ellsworth, D. S., Prentice,
  I. C., Barton, C. V., ... & Wingate, L. (2011). Reconciling the
  optimal and empirical approaches to modelling stomatal conductance.
  Global Change Biology, 17(6), 2134-2144.
* Wang, Y. P., Kowalczyk, E., Leuning, R., Abramowitz, G., Raupach,
  M. R., Pak, B., ... & Luhar, A. (2011). Diagnosing errors in a land
  surface model (CABLE) in the time and frequency domains. Journal of
  Geophysical Research: Biogeosciences, 116(G1).

"""

__title__ = "Typical LSM An & E iterative solving with the USO model"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (19.02.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import numpy as np  # array manipulations, math operators

# own modules
from TractLSM import conv, cst  # unit converter & general constants
from TractLSM.SPAC import fwsoil  # soil moisture stress
from TractLSM.SPAC import absorbed_radiation_2_leaves  # radiation
from TractLSM.SPAC import conductances, leaf_temperature  # energy
from TractLSM.SPAC import LH_water_vapour, vpsat, slope_vpsat, psychometric
from TractLSM.SPAC import calc_photosynthesis, rubisco_limit  # physio


# ======================================================================

def calc_trans(p, Tleaf, gs):

    """
    Calculates transpiration following Penman-Monteith at the leaf level
    accounting for effects of leaf temperature and feedback on
    evaporation.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Tleaf: float
        leaf temperature [degC]

    gs: float
        stomatal conductance [mol m-2 s-1]

    Returns:
    --------
    trans: float
        transpiration rate [mol m-2 s-1]

    real_zero: boolean
        True if the transpiration is really zero, False if Rnet is
        negative

    gw: float
        total leaf conductance to water vapour [mol m-2 s-1]

    gb: float
        boundary layer conductance to water vapour [mol m-2 s-1]

    """

    # check that the trans value satisfies the energy balance
    real_zero = True

    # get conductances
    gw, gH, gb, __ = conductances(p, Tleaf=Tleaf, gs=gs)  # mol m-2 s-1

    # latent heat of water vapor
    Lambda = LH_water_vapour(p)  # J mol-1

    # slope of saturation vapour pressure of water vs Tair
    slp = slope_vpsat(p)  # kPa degK-1

    if np.isclose(gs, 0., rtol=cst.zero, atol=cst.zero):
        trans = cst.zero

    else:
        gamm = psychometric(p)  # psychrometric constant, kPa degK-1
        trans = (slp * p.Rnet + p.VPD * gH * cst.Cp) / (Lambda *
                                                        (slp + gamm * gH / gw))

        if trans < 0.:  # Penman-Monteith failed, non-physical trans
            real_zero = False

        trans = max(cst.zero, trans)  # mol m-2 s-1

    return trans, real_zero, gw, gb


def solve_std(p, sw, photo='Farquhar', threshold_conv=0.015, iter_max=40):

    """
    Checks the energy balance by looking for convergence of the new leaf
    temperature with the leaf temperature predicted by the previous
    iteration. Then returns the corresponding An, E, Ci, etc.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    sw: float
        mean volumetric soil moisture content [m3 m-3]

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
    trans_can: float
        transpiration rate of canopy [mmol m-2 s-1] across leaves

    gs_can: float
        stomatal conductance of canopy [mol m-2 s-1] across leaves

    An_can: float
        C assimilation rate of canopy [umol m-2 s-1] across leaves

    Ci_can: float
        average intercellular CO2 concentration of canopy [Pa] across
        leaves

    rublim_can: string
        'True' if the C assimilation is rubisco limited, 'False'
        otherwise.

    """

    # retrieve sunlit / shaded fractions
    fRcan, fPPFD, fLAI, fscale2can, fgradis = absorbed_radiation_2_leaves(p)

    # saturation vapour pressure of water at Tair
    esat_a = vpsat(p.Tair)  # kPa

    # sunlit / shaded outputs
    trans = np.zeros(len(fPPFD))
    gs = np.zeros(len(fPPFD))
    An = np.zeros(len(fPPFD))
    Aj = np.zeros(len(fPPFD))
    Ac = np.zeros(len(fPPFD))
    Ci = np.zeros(len(fPPFD))

    # original LAI, PPFD
    LAI = p.LAI
    PPFD = p.PPFD

    # sunlit / shaded loop
    for i in range(len(fRcan)):

        p.Rnet = fRcan[i]
        p.PPFD = fPPFD[i]
        p.LAI = fLAI[i]
        p.scale2can = fscale2can[i]
        p.gradis = fgradis[i]

        if p.PPFD > 50.:  # min threshold for photosynthesis
            fw = fwsoil(p, sw)
            Cs = p.CO2  # Pa
            Tleaf = p.Tair  # deg C
            Dleaf = np.maximum(0.05, p.VPD)  # gs model not valid < 0.05

            # initialise gs_over_A
            g0 = 1.e-9  # g0 ~ 0, removing it entirely introduces errors
            Cs_umol_mol = Cs * conv.MILI * conv.FROM_kPa  # umol mol-1
            gs_over_A = g0 + (1. + p.g1 * fw / (Dleaf ** 0.5)) / Cs_umol_mol

            # iter on the solution until it is stable enough
            iter = 0

            while True:

                An[i], Aj[i], Ac[i], Ci[i] = \
                    calc_photosynthesis(p, 0., Cs, photo, Tleaf=Tleaf,
                                        gs_over_A=gs_over_A)

                # stomatal conductance, with fwsoil effect
                Cs_umol_mol = Cs * conv.MILI * conv.FROM_kPa
                gs_over_A = (g0 + (1. + p.g1 * fw / (Dleaf ** 0.5)) /
                             Cs_umol_mol)
                gs[i] = np.maximum(cst.zero, conv.GwvGc * gs_over_A * An[i])

                # calculate new trans, gw, gb, Tleaf
                trans[i], real_zero, gw, gb = calc_trans(p, Tleaf, gs[i])
                new_Tleaf, __ = leaf_temperature(p, trans[i], Tleaf=Tleaf,
                                                 gradis=True)

                # new Cs (in Pa)
                boundary_CO2 = (conv.ref_kPa * conv.FROM_MILI * An[i] /
                                (gb * conv.GbcvGb + gs[i] * conv.GcvGw))
                Cs = np.maximum(cst.zero,
                                np.minimum(p.CO2, p.CO2 - boundary_CO2))

                if (np.isclose(trans[i], cst.zero, rtol=cst.zero,
                   atol=cst.zero) or np.isclose(gw, cst.zero, rtol=cst.zero,
                   atol=cst.zero) or np.isclose(gs[i], cst.zero, rtol=cst.zero,
                   atol=cst.zero)):
                    Dleaf = np.maximum(0.05, p.VPD)  # kPa

                else:
                    esat_l = vpsat(new_Tleaf)  # vpsat at new Tleaf, kPa
                    Dleaf = (esat_l - (esat_a - p.VPD))  # leaf-air vpd, kPa

                # force stop when atm. conditions yield E < 0. (non-physical)
                if (iter < 1) and (not real_zero):
                    real_zero = None

                # check for convergence
                if ((real_zero is None) or (iter > iter_max) or ((real_zero)
                   and (abs(Tleaf - new_Tleaf) <= threshold_conv) and not
                   np.isclose(gs[i], cst.zero, rtol=cst.zero, atol=cst.zero))):
                    break

                # no convergence, iterate on leaf temperature
                Tleaf = new_Tleaf
                iter += 1

            if (np.isclose(trans[i], cst.zero, rtol=cst.zero, atol=cst.zero) or
               np.isclose(Ci[i], 0., rtol=cst.zero, atol=cst.zero) or
               (Ci[i] < 0.) or np.isclose(Ci[i], p.CO2, rtol=cst.zero,
                                          atol=cst.zero) or
               (Ci[i] > p.CO2) or (real_zero is None) or (not real_zero)):
                trans[i], gs[i], An[i], Ci[i] = (0.,) * 4

        else:
            trans[i], gs[i], An[i], Aj[i], Ac[i], Ci[i] = (0.,) * 6

    # scale to canopy: sum contributions from sunlit and shaded leaves
    with np.errstate(invalid='ignore'):  # if nans, do not raise warning
        trans_can = np.nansum(trans) * conv.MILI  # mmol m-2 s-1
        gs_can = np.nansum(gs)  # mol m-2 s-1
        An_can = np.nansum(An)  # umol m-2 s-1
        Ci_can = np.nanmean(Ci)  # Pa
        rublim_can = rubisco_limit(np.nansum(Aj), np.nansum(Ac))  # lim?

    # reset original LAI, PPFD
    p.LAI = LAI
    p.PPFD = PPFD

    if any(np.isnan([trans_can, gs_can, An_can, Ci_can])):
        trans_can, gs_can, An_can, Ci_can = (0.,) * 4

    return trans_can, gs_can, An_can, Ci_can, rublim_can


def set_trans_std(p, daily_steps, sw, photo='Farquhar'):

    """
    Forces the An, E, and Ci in cases where the next solving time step
    leads to wilting of the plant / depletion of the soil model.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    daily_steps: int
        number of timesteps in a day

    sw: float
        time step's mean volumetric soil moisture content [m3 m-3]

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    Returns:
    --------
    trans_can: float
        transpiration rate of canopy [mmol m-2 s-1] across leaves

    gs_can: float
        stomatal conductance of canopy [mol m-2 s-1] across leaves

    An_can: float
        C assimilation rate of canopy [umol m-2 s-1] across leaves

    Ci: float
        average intercellular CO2 concentration of canopy [Pa] across
        leaves

    rublim_can: string
        'True' if the C assimilation is rubisco limited, 'False'
        otherwise.

    """

    # unit conversions
    mm_2_molsqrtm = cst.rho / cst.MH2O  # mm to mol m-2

    # soil moisture conditions
    fw = fwsoil(p, sw)

    if np.isclose(fw, 0.):  # if fw too small, cannot solve for An
        return (0.,) * 5

    # leftover water that can be transpired by the plant (mol m-2 s-1)
    avail_H2O = (sw - p.pwp) * (p.Zbottom - p.Ztop) * conv.MILI
    trans_can = avail_H2O * mm_2_molsqrtm * conv.HLFHR_2_SEC

    # approximate Ci (Medlyn et al., 2011)
    Ci = np.maximum(cst.zero, p.CO2 * p.g1 * fw / (p.g1 * fw *
                    (np.maximum(0.05, p.VPD) ** 0.5)))

    # approximate An
    An_can, Aj_can, Ac_can = calc_photosynthesis(p, trans_can, Ci, photo,
                                                 Tleaf=p.Tair)

    # approximate gs (mol m-2 s-1)
    g0 = 1.e-9  # g0 ~ 0, but removing it entirely introduces errors
    Ci_umol_mol = Ci * conv.MILI * conv.FROM_kPa  # Pa to umol mol-1
    gs_over_A = (g0 + (1. + p.g1 * fw / (np.maximum(0.05, p.VPD) ** 0.5)) /
                 Ci_umol_mol)
    gs_can = np.maximum(cst.zero, conv.GwvGc * gs_over_A * An_can)

    if ((An_can <= cst.zero) or np.isclose(Ci, 0., rtol=cst.zero,
       atol=cst.zero) or (Ci < 0.) or np.isclose(Ci, p.CO2, rtol=cst.zero,
       atol=cst.zero) or (Ci > p.CO2)):  # cannot solve
        return (0.,) * 5

    rublim_can = rubisco_limit(Aj_can, Ac_can)  # which limitation?

    if any(np.isnan([trans_can, gs_can, An_can, Ci])):
        trans_can, gs_can, An_can, Ci = (0.,) * 4

    return trans_can * conv.MILI, gs_can, An_can, Ci, rublim_can
