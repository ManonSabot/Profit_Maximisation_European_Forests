# -*- coding: utf-8 -*-

"""
Get the solar zenith angle: this is used to calculate the canopy diffuse
and direct radiation in the two-leaf model.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Bonan, Gordon B. "A simulation model of environmental processes and
  vegetation patterns in Boreal forests: Test case Fairbanks, Alaska."
  (1988).
* Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
  biophysics. Pg 169.
* De Pury & Farquhar (1997) PCE, 20, 537-557.
* Hughes, David W.; Yallop, B. D.; Hohenkerk, C. Y. (1989), "The
  Equation of Time", Monthly Notices of the Royal Astronomical Society
  238: 1529–1535
* Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
* Nikolov, Nedialko T., and Karl F. Zeller. "A solar radiation algorithm
  for ecosystem dynamic models." Ecological Modelling 61.3-4 (1992):
  149-168.
* Spencer, J. W. (1971). Fourier series representation of the position
  of the sun.

"""

__title__ = "solar geometry of location"
__author__ = ["Manon E. B. Sabot", "Martin De Kauwe"]
__version__ = "2.0 (12.11.2018)"
__email__ = ["m.e.b.sabot@gmail.com", "mdekauwe@gmail.com"]


# ======================================================================

# general modules
import numpy as np  # array manipulations, math operators

# own modules
try:
    from constants_and_conversions import ConvertUnits
    from constants_and_conversions import Constants

except (ImportError, ModuleNotFoundError):
    from TractLSM.Utils.constants_and_conversions import ConvertUnits
    from TractLSM.Utils.constants_and_conversions import Constants


# ======================================================================

def day_angle(doy):

    """
    Calculation of day angle. De Pury & Farquhar eq. A18

    Returns:
    ---------
    fractional year / day angle [radians]

    """
    return 2. * np.pi * doy / 365.25


def solar_declination(doy):

    """
    Solar Declination Angle is a function of day of year and is
    indepenent of location, varying between 23deg45' and -23deg45'. De
    Pury & Farquhar eq. A14

    Arguments:
    ----------
    doy: int or float
        day of year

    gamma: float
        fractional year / day angle [radians]

    Returns:
    --------
    Solar Declination Angle [radians]

    """

    return np.radians(-23.45) * np.cos(2. * np.pi * (doy + 10.) / 365.25)


def eq_of_time(gamma):

    """
    Equation of time - correction for the difference btw solar time and
    the clock time.

    Arguments:
    ----------
    doy: int of float
        day of year

    gamma: float
        fractional year / day angle [radians]

    Returns:
    --------
    Equation of time [minutes]

    """

    # Spencer '71. This best matches the de Pury worked example (pg 554)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(gamma) -
                    0.032077 * np.sin(gamma) - 0.014615 * np.cos(2. * gamma) -
                    0.04089 * np.sin(2. * gamma))

    return eqt


def solar_midday(eqt, longitude):

    """
    Calculation of solar midday. De Pury & Farquhar eq. A16

    Arguments:
    ----------
    eqt: float
        equation of time [minutes]

    longitude: float
        longitude [deg]

    Returns:
    --------
    Solar midday [hours]

    """

    # international standard meridians are multiples of 15 deg E/W
    lonmed = round(longitude / 15.) * 15.

    return 12. + (4. * (lonmed - longitude) - eqt) / 60.


def hour_angle(hod, gamma, longitude):

    """
    Calculation of hour angle, relative to solar midday. De Pury &
    Farquhar eq. A15

    Arguments:
    ----------
    hod: float
        hour of the day (0.5 to 24)

    gamma: float
        fractional year / day angle [radians]

    longitude: float
        longitude [deg]

    Returns:
    ---------
    Hour angle [radians]

    """

    eqt = eq_of_time(gamma)
    t0 = solar_midday(eqt, longitude)

    return np.pi * (hod - t0) / 12.


def cos_zenith(doy, hod, latitude, longitude):

    """
    The solar zenith angle is the angle between the zenith and the
    centre of the sun's disc. The solar elevation angle is the altitude
    of the sun, the angle between the horizon and the centre of the
    sun's disc. Since these two angles are complementary, the cosine of
    either one of them equals the sine of the other, i.e.
    cos theta = sin beta.

    Arguments:
    ----------
    doy: int or float
        day of year

    hod: float
        hour of the day (0.5 to 24)

    latitude: float
        latitude [deg]

    longitude: float
        longitude [deg]

    Returns:
    --------
    cosine of the zenith angle of the sun (0-1)

    """

    gamma = day_angle(doy)
    rdec = solar_declination(doy)
    h = hour_angle(hod, gamma, longitude)

    # A13 - De Pury & Farquhar
    cos_zenith = (np.sin(np.radians(latitude)) * np.sin(rdec) +
                  np.cos(np.radians(latitude)) * np.cos(rdec) * np.cos(h))

    return np.maximum(0., np.minimum(1., cos_zenith))


def top_atm_rad(doy, cos_zen):

    """
    Solar radiation incident outside the earth's atmosphere, e.g.
    extra-terrestrial radiation. The value varies a little with the
    earth's orbit. Using formula (eq 1) from Spitters et al. (1986).

    Arguments:
    ----------
    doy: int or array
        day of year

    cos_zen: array
        cosine of the zenith angle of the sun (0-1)

    Returns:
    --------
    Sxt: float
        solar radiation incident at the top of the atmosphere
        [umol m-2 s-1]

    """

    Sxt = Constants().S0 * cos_zen * (1. + 0.033 * np.cos(2. * np.pi * doy /
                                                          365.25))

    try:
        Sxt[cos_zen <= 1.e-10] = 0.

    except TypeError:
        if cos_zen <= 1.e-10:
            Sxt = 0.

    return Sxt * ConvertUnits().SW_2_PAR  # umol m-2 s-1


def cloud_correction(Sxt, cld_cov, latitude):

    """
    We can account for scattering and absorbing effects of air
    molecules, water vapor, dust, ozone, and carbon dioxide as a linear
    function of undepleted solar radiation and the average portion of
    the sky covered by clouds, as in Bonan (1988) and Nikolov and
    Zeller (1992).

    Arguments:
    ----------
    Sxt: float
        solar radiation incident at the top of the atmosphere
        [umol m-2 s-1]

    cld_cov: float
        cloud cover (0-1)

    latitude: float
        latitude [deg]

    Returns:
    --------
    The solar radiation near the bottom of the boundary layer, corrected
    for clouds [umol m-2 s-1]

    """

    # empirical correction parameters; beta and sigma are unitless
    costerm = 1. - 1.3614 * np.cos(np.radians(latitude))
    alpha = 32.9835 - 64.884 * costerm  # in cal cm-2 d-1
    beta = 0.715 - 0.3183 * costerm
    sigma = 0.03259 + 0.02 * np.log(np.maximum(cld_cov * 10., 0.001))

    # to apply the empirical corrections, alpha must be in umol m-2 s-1
    alpha *= (ConvertUnits().Cal_2_J * ConvertUnits().J_2_UMOL *
              ConvertUnits().MILI * 10. / ConvertUnits().SEC_2_DAY)

    # empirical corrections
    Sxt *= beta - sigma * cld_cov
    Sxt -= alpha
    Sxt[Sxt < 0.] = 0.

    return Sxt


def spitters(cos_zenith, Sxt, PPFD):

    """
    Algorithm to estimate the diffuse component from the measured
    irradiance as given in eq 20a-d of Spitters et al., 1986.

    Arguments:
    ----------
    cos_zenith: array
        cosine of the zenith angle of the sun (0-1)

    Sxt: array
        solar radiation incident at the top of the atmosphere
        [umol m-2 s-1]

    PPFD: array
        PAR [umol m-2 s-1]

    Returns:
    --------
    direct_frac: array
        the direct component of incoming radiation

    diffuse_frac: array
        the diffuse component of incoming radiation

    """

    # atmospheric transmisivity
    tau = PPFD / Sxt
    tau[cos_zenith <= 1.e-10] = 0.
    tau[PPFD * ConvertUnits().PAR_2_SW <= 10.] = 0.

    # Spitters formula
    R = 0.847 + cos_zenith * (1.04 * cos_zenith - 1.61)

    # diffuse fraction
    diffuse_frac = R
    diffuse_frac[np.logical_and(tau > 0.22, tau <= 0.35)] = \
        (1. - 6.4 * (tau[np.logical_and(tau > 0.22, tau <= 0.35)] -
                     0.22) ** 2.)
    diffuse_frac[np.logical_and(tau > 0.35, tau <= (1.47 - R) / 1.66)] = \
        1.4728 - 1.66 * tau[np.logical_and(tau > 0.35, tau <= (1.47 - R) /
                                           1.66)]
    diffuse_frac[cos_zenith <= 1.e-2] = 1.
    diffuse_frac = np.maximum(0., np.minimum(1., diffuse_frac))

    return 1. - diffuse_frac, diffuse_frac


def LAI_sunlit(cos_zenith, direct_frac, LAI, LAI_thresh):

    """
    Calculate the sunlit LAI fraction of the canopy. The Ross-Goudriaan
    algorithm is used to get the ratio of the projected area of leaves
    in the direction perpendicular to the direction of incident solar
    radiation and the actual leaf area (Sellers (1985) eq. 13;
    CABLE code, Kowalczyk (2006) eq. 28/29).

    Arguments:
    ----------
    cos_zenith: array
        cosine of the zenith angle of the sun (0-1)

    direct_frac: array
        the direct component of incoming radiation

    LAI: float or array
        leaf area index [m2 m-2]

    LAI_thresh: bool or array
        mask for LAI values meeting the threshold [m2 m-2]

    Returns:
    --------
    LAI_sun: float or array
        LAI of the sunlit leaves [m2 m-2]

    """

    # set values for use in the function
    cos3 = np.array([np.cos(np.deg2rad(15.)), np.cos(np.deg2rad(45.)),
                     np.cos(np.deg2rad(75.))])  # 15, 45, 75 deg angles
    gauss_w = np.array([0.308, 0.514, 0.178])  # Gaussian integ. weights

    # leaf angle parameters
    chi_l = 9.99999978E-03  # leaf angle dist related (=0 for spherical)
    xphi1 = 0.5 - chi_l * (0.633 + 0.33 * chi_l)
    xphi2 = 0.877 * (1. - 2. * xphi1)

    # Ross-Goudriaan's algo (approx as eq 28; Kowalcyk et al., 2006)
    RG = xphi1 + xphi2 * cos_zenith

    # black leaves dir beam extinct coef. (eq 26; Kowalcyk et al., 2006)
    kb = RG / cos_zenith  # veg, light
    kb[np.logical_and(LAI_thresh, direct_frac <= 1.e-3)] = 0.5

    # black leaves diff extinct coef. (eq 27; Kowalcyk et al., 2006)
    kbx = (xphi1 + xphi2 * cos3) / cos3  # approximate kb integration
    integration_term = (np.nan_to_num(gauss_w[0] * np.exp(-kbx[0] * LAI)) +
                        np.nan_to_num(gauss_w[1] * np.exp(-kbx[1] * LAI)) +
                        np.nan_to_num(gauss_w[2] * np.exp(-kbx[2] * LAI)))
    kd = -np.log(integration_term) / LAI
    kd[LAI_thresh] = 0.7  # i.e. bare soil

    # radiation thresholds
    kb[np.abs(kb - kd) < 1.e-3] = kd[np.abs(kb - kd) < 1.e-3] + 1.e-3
    kb[direct_frac < 1.e-3] = 1.e5

    # sunlit & shaded LAI frac of can (eq 18; De Pury & Farquhar, 1997)
    LAI_sun = (1. - np.exp(-kb * LAI)) / kb

    return LAI_sun


def composite_LAI(doy, cos_zenith, Sxt, PPFD, LAI):

    """
    Scatter the LAI into the sun/shade fractions and then weight it to
    get a reasonable composite LAI.

    Arguments:
    ----------
    doy: array
        day of year

    cos_zenith: array
        cosine of the zenith angle of the sun (0-1)

    Sxt: array
        solar radiation incident at the top of the atmosphere
        [umol m-2 s-1]

    PPFD: array
        PAR [umol m-2 s-1]

    LAI: float or array
        leaf area index [m2 m-2]

    Returns:
    --------
    LAI: float
        composite LAI [m2 m-2]

    """

    # direct/diffuse light fractions
    fdir, fdiff = spitters(cos_zenith, Sxt, PPFD)

    # valid thresholds
    LAI_thresh = (LAI <= Constants().LAI_thresh)
    RAD_thresh = Constants().RAD_thresh

    # LAI fractions
    LAI_sun = LAI_sunlit(cos_zenith, fdir, LAI, LAI_thresh)
    LAI_sun[np.logical_and(LAI_thresh,
                           PPFD * ConvertUnits().PAR_2_SW <= RAD_thresh)] = 0.
    LAI_sha = LAI - LAI_sun

    # new composite LAI, keeping low PAR threshold
    LAI_comp = np.ma.masked_where(LAI_sun <= 0.,
                                  fdiff * LAI_sha + fdir * LAI_sun)
    return LAI_comp
