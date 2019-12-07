#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Default parameter class, necessary to run the model.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* ten Berge, Hein FM. Heat and water transfer in bare topsoil and lower
  atmosphere. No. 33. Pudoc, 1990.
* Campbell, G. S., & Norman, J. M. “An Introduction to Environmental
  Biophysics” 2nd Edition, Springer-Verlag, New York, 1998.
* Choat, B., Jansen, S., Brodribb, T. J., Cochard, H., Delzon, S.,
  Bhaskar, R., ... & Jacobsen, A. L. (2012). Global convergence in the
  vulnerability of forests to drought. Nature, 491(7426), 752.
* Kattge, J., & Knorr, W. (2007). Temperature acclimation in a
  biochemical model of photosynthesis: a reanalysis of data from 36
  species. Plant, cell & environment, 30(9), 1176-1190.
* Medlyn, B. E., Dreyer, E., Ellsworth, D., Forstreuter, M., Harley,
  P. C., Kirschbaum, M. U. F., ... & Wang, K. (2002). Temperature
  response of parameters of a biochemically based model of
  photosynthesis. II. A review of experimental data. Plant, Cell &
  Environment, 25(9), 1167-1179.
* Peltoniemi, M. S., Duursma, R. A., & Medlyn, B. E. (2012). Co-optimal
  distribution of leaf nitrogen and hydraulic conductance in plant
  canopies. Tree Physiology, 32(5), 510-519.

"""

__title__ = "default parameter class necessary to run the model"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (24.01.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

class default_params(object):  # default inputs needed to run model

    def __init__(self):

        # information used by the weather generator
        self.doy = 180.  # day of the year
        self.tmin = 2.  # degC
        self.tmax = 24.
        self.rain_day = 2.  # mm d-1
        self.vpd15prev = 3.4
        self.vpd09 = 1.4
        self.vpd15 = 2.3
        self.vpd09next = 1.8
        self.sw_rad_day = 1080. * 10.  # W m-2, for 10 daylight hours
        self.Patm = 101.325  # kPa
        self.u = 2.  # m s-1

        # location is generally important!
        self.lat = 38.569120
        self.lon = -80.018519

        # gas concentrations
        self.CO2 = 37.  # Pa, by default 375 ppm as in ~2005
        self.O2 = 20.73  # kPa

        # photosynthesis related
        self.Vmax25 = 100.  # max carboxyl rate @ 25 degC (umol m-2 s-1)
        self.gamstar25 = 4.22  # CO2 compensation point @ 25 degC (Pa)
        self.Tref = 25.  # ref T for Vmax25, gamstar, deltaS, Hd
        self.JV = 1.67  # Jmax25 to Vmax25 ratio (Medlyn et al., 2002)
        self.Rlref = self.Vmax25 * 0.015  # resp @ TRlref (umol m-2 s-1)
        self.TRlref = 25.  # T for the ref keaf respiration (degC)
        self.Kc25 = 39.96  # Michaelis-Menten cst for carboxylation (Pa)
        self.Ko25 = 27.48  # Michaelis-Menten cst for oxygenation (kPa)
        self.alpha = 0.3  # quantum yield photo (mol(photon) mol(e-)-1)

        # Farquhar model
        self.c1 = 0.7  # curvature of light response
        self.c2 = 0.99  # transition Je vs Jc (Peltoniemi et al., 2012)

        # Collatz model
        self.c3 = 0.83  # curvature of light response
        self.c4 = 0.93  # transition Je vs Jc

        # energies of activation
        self.Ev = 60000.  # Vcmax, J mol-1
        self.Ej = 30000.  # Jmax, J mol-1
        self.Egamstar = 37830.  # gamstar, J mol-1
        self.Ec = 79430.  # carboxylation, J mol-1
        self.Eo = 36380.  # oxygenation, J mol-1

        # inhibition at higher temperatures (Kattge & Knorr, 2007)
        self.deltaSv = 650.  # Vmax entropy factor (J mol-1 K-1)
        self.deltaSj = 650.  # Jmax entropy factor (J mol-1 K-1)
        self.Hdv = 200000.  # Vmax decrease rate above opt T (J mol-1)
        self.Hdj = 200000.  # Jmax decrease rate above opt T (J mol-1)

        # relating to light / rad (C & N is Campbell & Norman)
        self.eps_l = 0.97  # leaf emiss. LW (Table 11.3 C & N 1998)
        self.eps_s = 0.945  # soil emiss. LW (Table 11.3 C & N 1998)
        self.albedo_l = 0.062  # leaf SW vis (CABLE)
        self.albedo_ws = 0.1  # wet soil SW vis (Table 11.2 C & N 1998)
        self.albedo_ds = 0.25  # dry soil SW vis (ten Berge, 1990)
        self.tau_l = 0.05  # leaf transmis. SW vis (CABLE)
        self.chi_l = 9.99999978E-03  # leaf angle dist (spherical = 0)
        self.kn = 0.001  # extinction coef. of nitrogren (CABLE)

        # hydraulics / plant
        self.LAI = 2.  # m2 m-2
        self.max_leaf_width = 0.001  # m
        self.canopy_intercept = 0.05  # % rainfall intercepted by canopy
        self.P50 = 6.6  # xylem pressure at P50 (-MPa) - J. virginiana
        self.P88 = 10.5  # same at P88 (-MPa) (Choat et al., 2012)
        self.kmax = 1.  # max plant hydr cond / LAI (mmol m-2 s-1 MPa-1)
        self.ratiocrit = 0.05  # degree stom control? kcrit = N%(kmax)
        self.g1 = 2.35  # sensitivity of stomatal conduc to An (kPa0.5)
        self.sfw = 1.  # sensitivity power factor on fwsoil (unitless)
        self.nfc = 0.  # nudging term applied to the fc (m3 m-3)
        self.npwp = 0.  # nudging term applied to the pwp (m3 m-3)

        # soil
        self.ground_area = 1.  # m2
        self.Ztop = 0.02  # top soil layer depth (m)
        self.Zbottom = 1.  # bottom soil layer depth (m)
        self.Psie = -0.8e-3  # air entry point water potential (MPa)
        self.Ps = self.Psie  # initial soil water potential (MPa)
        self.hyds = 8.05e-6  # saturation hydraulic conductivity (m s-1)
        self.fc = 0.3  # field capacity (m3 m-3)
        self.theta_sat = 0.5  # soil moisture at saturation (m3 m-3)
        self.pwp = 0.05  # plant wilting point (m3 m-3)
        self.bch = 4.  # Clapp and Hornberger index
        self.r_soil = 0.6  # soil resistance to water infil (unitless)

        return
