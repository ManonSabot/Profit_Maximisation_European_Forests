# -*- coding: utf-8 -*-

"""
Simple tipping bucket soil water balance model, a proxy for soil
hydrology.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

Reference:
-----------
* Clapp, R. B., & Hornberger, G. M. (1978). Empirical equations for some
  soil hydraulic properties. Water resources research, 14(4), 601-604.
* Cosby, B. J., Hornberger, G. M., Clapp, R. B., & Ginn, T. (1984). A
  statistical exploration of the relationships of soil moisture
  characteristics to the physical properties of soils. Water resources
  research, 20(6), 682-690.
* Duursma et al. (2008). Predicting the decline in daily maximum
  transpiration rate of two pine stands during drought based on constant
  minimum leaf water potential and plant hydraulic conductance. Tree
  Physiology, 28(2), 265-276.
* Monteith, J. L., & Unsworth, M. H. (1990). Principles of environmental
  physics. Arnold. SE, London, UK.
* Ritchie, J. T. (1972). Model for predicting evaporation from a row
  crop with incomplete cover. Water resources research, 8(5), 1204-1213.

"""

__title__ = "Tipping bucket soil water module"
__author__ = "Manon E. B. Sabot"
__version__ = "2.0 (13.02.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import numpy as np  # array manipulations, math operators

# own modules
from TractLSM import conv, cst  # unit converter & general constants
from TractLSM.SPAC import slope_vpsat, LH_water_vapour, psychometric
from TractLSM.SPAC import net_radiation  # canopy radiation


# ======================================================================

def fwsoil(p, sw):

    """
    Calculates the empirical soil moisture stress factor that determines
    the stomatal conductance and soil's responses to water limitation.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    sw: float
        mean volumetric soil moisture content [m3 m-3]

    Returns:
    --------
    The empirical stomatal conductance and soil's responses to soil
    moisture stress.

    """

    fw = (np.maximum(cst.zero, sw - (p.pwp + p.npwp)) /
          np.maximum(cst.zero, (p.fc + p.nfc) - (p.pwp + p.npwp))) ** p.sfw

    return np.maximum(cst.zero, np.minimum(1., fw))


def evap_bare_soil(p):

    """
    Evaporation at the potential/equilibrium rate, where aerodynamic
    conductance is zero (i.e. winds are calm).

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    Returns:
    --------
    The evaporation of bare soil [mol m-2 s-1], using Penman eq of
    Monteith and Unsworth, 1990.

    """

    # latent heat of water vapor
    Lambda = LH_water_vapour(p)  # J mol-1

    # psychrometric constant
    gamm = psychometric(p)  # kPa deg-1

    # slope of saturation vapour pressure of water vs Tair
    slp = slope_vpsat(p)  # kPa deg-1

    # net radiation of a surface
    Rnet = net_radiation(p, surface='soil')  # W m-2

    return np.maximum(0., (1. - p.r_soil) * Rnet * slp / ((slp + gamm) *
                                                          Lambda))


def soil_evap(p, sw):

    """
    Use Penman eq to calculate top soil evaporation flux. Soil
    evaporation is dependent upon soil wetness and plant cover. The net
    radiation term is scaled for the canopy cover and the impact of soil
    wetness is accounted for. As the soil dries the evaporation
    component reduces significantly.

    Key assumptions from Ritchie, 1972...
    When plant provides shade for the soil surface, evaporation will not
    be the same as bare soil evaporation. Wind speed, net radiation and
    VPD will all be lowered in proportion to the canopy density.
    Following Ritchie role of wind, VPD are assumed to be negligible and
    are therefore ignored. These assumptions are based on work with
    crops, where a fit is formed between the LAI of 5 crops types and
    the fraction of observed net radiation at the surface. Whilst the
    LAI does cover a large range, nominal 0–6, there are only 12
    measurements and only three from LAI > 3.and whether this holds for
    tree shading where the height from the soil to the base of the crown
    is larger is questionable.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    sw: float
        volumetric soil water content available for evaporation [m3 m-3]

    Returns:
    --------
    The evaporation of covered soil [mmol m-2 s-1], using Ritchie (1972)
    empirical LAI fit.

    """

    # account for LAI cover
    evap = evap_bare_soil(p) * np.exp(-0.398 * p.LAI)

    # account for soil wetness state
    evap *= fwsoil(p, sw) * conv.MILI  # mmol m-2 s-1

    if (np.isclose(evap, 0.)) or (evap < 0.):
        evap = 0.

    return evap


def drainage(p, daily_steps, sw, volume, depth):

    """
    Updates water content in the top and lower layer of the bucket by
    applying a drainage rate, with the low layer being split in 5
    different sublayers.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    daily_steps: int
        number of timesteps in a day

    sw: array
        volume of soil water in each layer [m3]

    volume: array
        volume of each layer [m3]

    depth: array
        depth of each layer [m]

    Returns:
    --------
    sw: array
        volume of soil water in each layer [m3]

    """

    # volumetric soil water contents
    sw /= volume  # m3 m-3

    # difference in soil water between consecutive layers
    delta = np.ediff1d(sw)  # m3 m-3

    # drainage speed, last layer is limited by saturation
    speed = (p.hyds * (np.minimum(sw[:-1], sw[1:]) / p.theta_sat)
             ** (2. * p.bch + 2.))  # m s-1
    speed = np.append(speed, p.hyds * (np.minimum(sw[-1], p.theta_sat) /
                                       p.theta_sat) ** (2. * p.bch + 2.))

    # limit speed for stability, m s-1
    time_spent = (conv.SEC_2_DAY / daily_steps)
    speed[:-1] = np.minimum(speed[:-1], 0.5 * depth[:-1] / time_spent)
    speed[-1] = np.minimum(0.5 * speed[-1], 0.5 * depth[-1] / time_spent)

    # calculate the drainage flux for the upper layers
    ratio = delta[1:] / (delta[:-1] + cst.zero * np.sign(delta[:-1]))
    ratio = np.append(np.nan_to_num(ratio), 0.)  # m3 m-3
    phi = np.maximum(0., np.minimum(1., 2. * ratio))
    phi = np.maximum(phi, np.minimum(2., ratio))
    flux = speed[:-1] * (sw[:-1] + phi * (np.minimum(sw[:-1], sw[1:]) -
                         sw[:-1]))  # m s-1

    # drainage flux for the lowest layer
    flux = np.append(flux, np.maximum(0., speed[-1] * sw[-1]))  # m s-1

    # constrain fluxes by the saturation point
    flux[:-1] = np.minimum(flux[:-1], flux[1:] + (p.theta_sat - sw[1:]) *
                           depth[1:] / time_spent)

    # update water balance, constrain by saturation point
    sw[1:] = np.minimum(p.theta_sat, sw[1:] + time_spent * (flux[:-1] -
                        flux[1:]) / depth[1:])
    sw[0] -= time_spent * flux[0] / depth[0]

    # back to volumes of water
    sw *= volume  # m3

    return np.minimum(p.theta_sat * volume, sw)


def sublayer_runoff(sw, max_sw, in_sw):

    """
    Updates water content in the bucket by tipping incoming rainfall
    throughout the layers, which is a coarse estimate of below ground
    runoff.

    Arguments:
    ----------
    sw: array
        volume of soil water in each layer [m3]

    max_sw: array
        maximum volume of water in each layer (i.e. volume at
        saturation) [m3]

    in_sw: float
        incoming volume of water entering the layers [m3]

    Returns:
    --------
    sw: array
        volume of soil water in each layer [m3]

    """

    # infiltrate the sub layers with rain water
    for i in range(len(sw)):

        if in_sw > (max_sw[i] - sw[i]):  # the layer overfills!
            in_sw -= max_sw[i] - sw[i]  # less water enters next layer
            sw[i] = max_sw[i]  # saturated layer, m3

            if in_sw < 0.:  # no more infiltration

                return sw

        else:  # not all layers are saturated
            sw[i] += in_sw  # include remainder entering water, m3

    return sw


def wetness(p, daily_steps, sw0, sw1, sw2, sw3, sw4, sw5, Es, E, Tsoil):

    """
    Updates the simple bucket soil water balance by calculating the soil
    volumetric water content [m3 m-3].

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    daily_steps: int
        number of timesteps in a day

    sw0: float
        top soil layer volumetric soil moisture content from the
        previous timestep [m3 m-3]

    sw1: float
        #1 soil layer volumetric soil moisture content from the previous
        timestep [m3 m-3]

    sw2: float
        #2 soil layer volumetric soil moisture content from the previous
        timestep [m3 m-3]

    sw3: float
        #3 soil layer volumetric soil moisture content from the previous
        timestep [m3 m-3]

    sw4: float
        #4 soil layer volumetric soil moisture content from the previous
        timestep [m3 m-3]

    sw5: float
        #5 soil layer volumetric soil moisture content from the previous
        timestep [m3 m-3]

    Es: float
        soil evaporation rate from the previous timestep [mmol m-2 s-1]

    E: float
        canopy transpiration rate from the previous timestep
        [mmol m-2 s-1]

    Tsoil: float
        mean soil temperature

    Returns:
    --------
    sw: float
        volumetric soil water content [m3 m-3]

    sw0: float
        top soil layer volumetric soil moisture content [m3 m-3]

    sw1: float
        #1 soil layer volumetric soil moisture content [m3 m-3]

    sw2: float
        #2 soil layer volumetric soil moisture content [m3 m-3]

    sw3: float
        #3 soil layer volumetric soil moisture content [m3 m-3]

    sw4: float
        #4 soil layer volumetric soil moisture content [m3 m-3]

    sw5: float
        #5 soil layer volumetric soil moisture content [m3 m-3]

    sevap: float
        next step's soil evaporative rate [mmol m-2 s-1]

    """

    # unit conversions (to mm 1/2h-1)
    mmolsqrtm_2_mm = conv.FROM_MILI * cst.MH2O / cst.rho  # m3 mol-1
    Es *= mmolsqrtm_2_mm * conv.SEC_2_HLFHR
    E *= mmolsqrtm_2_mm * conv.SEC_2_HLFHR
    precip = (1. - p.canopy_intercept) * p.precip / conv.HLFHR_2_DAY

    # deal with non-half-hourly steps!
    if daily_steps != 48:
        Es *= 48. / daily_steps
        E *= 48. / daily_steps
        precip *= 48. / daily_steps

    # account for snow / frozen soil
    if p.Tair < cst.zero:
        precip = 0.

    # 5 sublayers, with relative depths matching CABLE % root access
    zsl = np.asarray([0.2885118652, 0.42331875403, 0.26203702018, 0.0260879855,
                      0.00004437506])  # normalised % root water access
    vol = np.insert((p.soil_volume - p.soil_top_volume) * zsl, 0,
                    p.soil_top_volume)  # m3
    depth = np.insert((p.Zbottom - p.Ztop) * zsl, 0, p.Ztop)  # m

    # volumes of soil water
    sw = np.asarray([sw0, sw1, sw2, sw3, sw4, sw5])  # m3 m-3
    sw *= vol  # m3

    # first, evaporate from the previous timestep
    sw[0] -= Es * conv.FROM_MILI * p.ground_area  # m3
    sw[0] = np.maximum(0., np.minimum(p.theta_sat * vol[0], sw[0]))

    # second, preferably transpire from the topmost layers
    left_over = E * conv.FROM_MILI * p.ground_area  # m3

    for i in range(len(sw)):  # loop over the layers

        if sw[i] > cst.zero:  # water is available
            sw[i] -= left_over  # m3, remove water from that layer

            if sw[i] < 0.:  # layer is empty, rm below
                left_over = -sw[i]  # m3

            else:  # do not rm below
                break

    # restrain to the valid range
    sw = np.maximum(0., np.minimum(p.theta_sat * vol, sw))

    if precip > 0.:  # third, deal with infiltrating rainfall
        sw = sublayer_runoff(sw, p.theta_sat * vol,
                             precip * conv.FROM_MILI * p.ground_area)
        sw = np.maximum(0., np.minimum(p.theta_sat * vol, sw))

    # fourth, apply drainage rate
    sw = drainage(p, daily_steps, sw, vol, depth)  # m3
    sw = np.maximum(0., np.minimum(p.theta_sat * vol, sw))

    # calculate the soil evap to be removed at beginning of next water balance
    if sw[0] > cst.zero:
        sevap = soil_evap(p, sw[0] / vol[0])  # mmol m-2 s-1

    else:
        sevap = 0.

    # account for snow / frozen soil
    if Tsoil < cst.zero:
        sevap = 0.

    # is it too much evap for the amount of soil water left?
    next_Es = (sevap * mmolsqrtm_2_mm * conv.SEC_2_HLFHR * conv.FROM_MILI *
               p.ground_area)  # m3

    if daily_steps != 48:
        next_Es *= 48. / daily_steps

    # if too much evap, adjust evap to all the water available
    if (sw[0] - next_Es) < 0.:
        sevap = (sw[0] * conv.HLFHR_2_SEC * conv.MILI /
                 (p.ground_area * mmolsqrtm_2_mm))  # mmol m-2 s-1

        # deal with non-half-hourly steps!
        if daily_steps != 48:
            sevap *= 48. / daily_steps

    # overall soil moisture content
    sw_all = np.sum(sw) / p.soil_volume  # m3 m-3

    # volumetric soil moisture contents
    sw /= vol  # m3 m-3

    return sw_all, sw[0], sw[1], sw[2], sw[3], sw[4], sw[5], sevap


def water_potential(p, sw):

    """
    Calculates the soil water potential [MPa]. The parameters bch and
    Psie are estimated using the Cosby et al. (1984) regression
    relationships from the soil sand/silt/clay fractions.

    Arguments:
    ----------
    p: pandas series
        time step's met data & params

    sw: float
        volumetric soil water content [m3 m-3]

    Returns:
    --------
    The soil water potential [MPa], Ps, using Clapp and Hornberger
    eq (1978)

    """

    if (sw is not None) and (sw >= cst.zero):
        return p.Psie * (sw / p.theta_sat) ** (- p.bch)

    elif np.isclose(abs(p.Ps), 0., rtol=cst.zero, atol=cst.zero):
        return p.Psie

    else:
        return p.theta_sat * (p.Ps / p.Psie) ** (-1. / p.bch)
