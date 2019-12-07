# -*- coding: utf-8 -*-

"""
Two classes for constants & regular unit conversions used throughout the
model's repository.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Asner, G. P., Scurlock, M. O., & Hicke, J. A. (2003). Global Ecology &
  Biogeography, 12, 191–205.
* McCree, K. J. (1972). Test of current definitions of
  photosynthetically active radiation against leaf photosynthesis data.
  Agricultural Meteorology, 10, 443-453.
* Annex 1 of http://www.fao.org/3/x0490e/x0490e00.htm

"""

__title__ = "Constants & regular unit conversions"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (24.01.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

class ConvertUnits(object):

    """
    Unit conversions for PAR, from J to umol, degC to degK... Also has
    ratios of unit water to unit heat or unit carbon for plant
    conductances, and temporal and magnitudinal conversions.

    """

    def __init__(self):

        # convs
        self.J_2_UMOL = 4.57  # J to umolquanta (400-700nm), McCree 1972
        self.UMOL_2_J = 1. / self.J_2_UMOL
        self.SW_2_PAR = 4.57 * 0.5  # SW (W m-2) to PAR (umol m-2 s-1)
        self.PAR_2_SW = 1. / self.SW_2_PAR
        self.Cal_2_J = 4.1868  # cal to J
        self.C_2_K = 273.15  # degC to degK
        self.ref_kPa = 101.325  # to ref kPa
        self.FROM_kPa = 1. / self.ref_kPa  # from kPa to anything

        # plant conductances
        self.GbhvGb = 0.93  # ratio of Gbheat:Gbwater
        self.GbvGbh = 1. / self.GbhvGb  # Gbwater:Gbheat
        self.GbhvGbc = 1.32  # Gbheat:GbCO2
        self.GbcvGbh = 1. / self.GbhvGbc  # GbCO2:Gbheat
        self.GbvGbc = self.GbvGbh * self.GbhvGbc  # Gbwater:GbCO2
        self.GbcvGb = 1. / self.GbvGbc  # GbCO2:Gbwater
        self.GwvGc = 1.57  # Gwater:GCO2
        self.GcvGw = 1. / self.GwvGc  # GCO2:Gwater

        # time
        self.HR_2_DAY = 24.
        self.SEC_2_HR = 60. * 60.
        self.SEC_2_DAY = self.SEC_2_HR * self.HR_2_DAY
        self.DAY_2_HR = 1. / self.HR_2_DAY
        self.DAY_2_SEC = 1. / self.SEC_2_DAY
        self.SEC_2_HLFHR = self.SEC_2_HR / 2.
        self.HLFHR_2_SEC = 1. / self.SEC_2_HLFHR
        self.HLFHR_2_DAY = self.SEC_2_DAY / self.SEC_2_HLFHR

        # magnitudes
        self.U = 1.e6
        self.FROM_U = 1.e-6
        self.MILI = 1.e3
        self.FROM_MILI = 1.e-3
        self.MEGA = 1.e-6
        self.FROM_MEGA = 1.e6

        # FAO equivalent water
        self.Wpm2_2_MJpm2pd = 0.0864
        self.MJpm2pd_2_mmpd = 0.408
        self.Wpm2_2_mmpday = self.Wpm2_2_MJpm2pd * self.MJpm2pd_2_mmpd
        self.Wpm2_2_mmphr = self.Wpm2_2_mmpday / self.HR_2_DAY
        self.Wpm2_2_mmphlfhr = self.Wpm2_2_mmphr / 2.

        # carbon conversions
        self.mmolH2Opm2ps_2_mmphr = (Constants().MH2O * self.FROM_U
                                                      * self.SEC_2_HR)
        self.mmolH2Opm2ps_2_mmphlfhr = self.mmolH2Opm2ps_2_mmphr / 2.
        self.umolCpm2ps_2_gCpm2phr = (Constants().MC * self.FROM_U
                                                     * self.SEC_2_HR)
        self.umolCpm2ps_2_gCpm2phlfhr = self.umolCpm2ps_2_gCpm2phr / 2.

        return


class Constants(object):

    """
    Constants for water and air properties, radiation, and leaf to
    anopy.

    """

    def __init__(self):

        self.zero = 1.e-17  # precision for the calculations, 17f

        self.g0 = 9.80665  # gravitational acceleration (m s-2)
        self.S0 = 1362.  # Solar constant (W m-2)
        self.vonK = 0.41  # Von Karman constant (-)
        self.Mair = 28.9644  # molar mass of dry air (g mol-1)
        self.MH = 1.00794  # molar mass of hydrogen (g mol-1)
        self.MC = 12.0107  # molar mass of C (g mol-1)
        self.MO = 15.999  # molar mass of O (g mol-1)
        self.MCO2 = self.MC + 2. * self.MO  # CO2 molar mass (g mol-1)
        self.MH2O = 2. * self.MH + self.MO  # H2O molar mass (g mol-1)
        self.rho = 1000.  # density of water (kg m-3)
        self.LH2O = 2.501e6  # latent heat H2O (J kg-1)
        self.DH = 21.5e-6  # molecular diffusivity to heat (m2 s-1)
        self.R = 8.3144598  # molar ideal gas constant (J mol-1 K-1)
        self.Rd = 286.9  # gas constant for dry air (J kg-1 K-1)
        self.Rv = 461.5  # gas constant for water vapor (J kg-1 K-1)
        self.Cp = 29.29  # heat capacity of dry air, cst P (J mol-1 K-1)
        self.Lb = 0.0065  # standard temperature lapse rate (K m-1)
        self.sigma = 5.67e-8  # Stefan-Boltzmann cst (W m-2 K-4)

        # threshold based on Asner et al., 2003
        self.LAI_thresh = 1.e-3  # minimum valid LAI (m2 m-2)
        self.RAD_thresh = 1.e-3  # minimum valid radiation (W m-2)

        return
