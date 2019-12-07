# -*- coding: utf-8 -*-

"""
Simple weather generator functions

References:
-----------
* De Pury, D. G. G., & Farquhar, G. D. (1997). Simple scaling of
  photosynthesis from leaves to canopies without the errors of big‐leaf
  models. Plant, Cell & Environment, 20(5), 537-557.
* Haverd, V., Raupach, M. R., Briggs, P. R., Canadell, J. G., Isaac, P.,
  Pickett-Heaps, C., ... & Wang, Z. (2013). Multiple observation types
  reduce uncertainty in Australia's terrestrial carbon and water cycles.
  Biogeosciences, 10(3), 2011-2040.
* Kimball, B. A., & Bellamy, L. A. (1986). Generation of diurnal solar
  radiation, temperature, and humidity patterns. Energy in agriculture,
  5(3), 185-197.
* Leuning, R., Kelliher, F. M., De Pury, D. G. G., & Schulze, E. D.
  (1995). Leaf nitrogen, photosynthesis, conductance and transpiration:
  scaling from leaves to canopies. Plant, Cell & Environment, 18(10),
  1183-1200.
* Loustau, D., Pluviaud, F., Bosc, A., Porté, A., Berbigier, P., Déqué,
  M., & Pérarnaud, V. (2001). Impact of a regional 2× CO2 climate
  scenario on the water balance, carbon balance and primary production
  of maritime pine in southwestern France. Models for the Sustainable
  Management of Plantation Forests. Ed. M. Tomé. European Cultivated
  Forest Inst., EFI Proc, (41D), 45-58.
* Parton, W. J., & Logan, J. A. (1981). A model for diurnal variation
  in soil and air temperature. Agricultural Meteorology, 23, 205-216.
* Spencer, J. W. (1971). Fourier series representation of the position
  of the sun. Search, 2(5), 172-172.
* Spitters, C. J. T., Toussaint, H. A. J. M., & Goudriaan, J. (1986).
  Separating the diffuse and direct component of global radiation and
  its implications for modeling canopy photosynthesis Part I. Components
  of incoming radiation. Agricultural and Forest Meteorology, 38(1-3),
  217-229.

"""

__title__ = "weather generator adapted from mdekauwe's AWAP version"
__reference__ = "https://github.com/mdekauwe/weather_generator"
__author__ = ["Manon E. B. Sabot", "Martin G. De Kauwe"]
__version__ = "1.0 (31.03.2017)"
__email__ = ["m.sabot@unsw.edu.au", "mdekauwe@gmail.com"]


# ======================================================================

# import general modules
import numpy as np  # array manipulations, math operators
from math import pi, cos, sin, exp, acos, asin  # math specs


# ======================================================================

def main(lat, lon, doy, tmin, tmax, rain_day, vpd15prev, vpd09, vpd15,
         vpd09next, sw_rad_day):

    """
    Main function: generates diurnal time courses of PAR, temperature,
                   precip, and vpd.

    Arguments:
    ----------
    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    doy: int
        day of year

    tmin: float
        day minimum temp [deg C]

    tmax: float
        day maximum temp [deg C]

    rain_day: float
        daily rainfall total [mm]

    vpd15_prev: float
        vpd at 3 pm the previous day [kPa]

    vpd09: float
        expected vpd at 9 am [kPa]

    vpd15: float
        expected vpd at 3 pm [kPa]

    vpd09_next: float
        expected vpd at 9 am the next day [kPa]

    sw_rad_day: float
        daily total short wave radiation [W m-2]

    Returns:
    ---------
    par: array
        diurnal course of the par [umol m-2 s-1]

    tday: array
        diurnal time course of the temperature [deg C]

    precip: array
        diurnal time course of rainfall [mm]

    vpd: array
        diurnal course of the vapor pressure deficit [kPa]

    """

    WG = WeatherGenerator(lat, lon)

    par_day = sw_rad_day * WG.SW_2_PAR
    par = WG.estimate_diurnal_par(par_day, doy)

    tday = WG.estimate_diurnal_temp(doy, tmin, tmax)

    precip = WG.disaggregate_rainfall(rain_day)

    vpd = WG.estimate_diurnal_vpd(vpd15prev, vpd09, vpd15, vpd09next)

    return par, tday, precip, vpd


# ======================================================================

# ~~~ Other functions are defined here ~~~

class WeatherGenerator(object):

    def __init__(self, lat, lon):

        """
        Unit conversions and initialisation of lat, lon, and the number
        of time steps in a day.

        Arguments:
        ----------
        lat: float
            latitude [decimal degrees]

        lon: float
            longitude [decimal degrees]

        """

        self.SW_2_PAR = 4.57 * 0.5  # SW (W m-2) to PAR (umol m-2 d-1)
        self.PAR_2_SW = 1. / self.SW_2_PAR
        self.SW_2_PAR_MJ = 0.5  # SW (MJ m-2 d-1) to PAR (umol m-2 d-1)
        self.J_TO_MJ = 1E-6
        self.SEC_2_DAY = 86400.
        self.MJ_TO_J = 1E6
        self.DAY_2_SEC = 1. / self.SEC_2_DAY
        self.SEC_2_HLFHR = 1800.
        self.HLFHR_2_SEC = 1. / self.SEC_2_HLFHR
        self.J_TO_UMOL = 4.57
        self.UMOL_TO_J = 1. / self.J_TO_UMOL
        self.UMOLPERJ = 4.57  # Conversion from J to umol quanta

        self.lat = lat
        self.lon = lon

        self.hours = np.arange(0., 24., 0.5)  # half hourly steps in day

        return

    def day_angle(self, doy):

        """
        Calculation of day angle - De Pury & Farquhar (1997): eq A18,
        also see Spencer (1971)

        Arguments:
        ----------
        doy: int
            day of year

        Returns:
        ---------
        the day angle [radians]

        """

        return 2. * pi * (float(doy) - 1.) / 365.

    def calculate_solar_declination(self, doy, gamma):

        """
        The Solar Declination Angle is a function of day of year and is
        indepenent of location. It varies between 23deg45' to -23deg45'.

        Arguments:
        ----------
        doy: int
            day of year

        gamma: float
            fractional year [radians]

        Returns:
        --------
        decl: float
            Solar Declination Angle [radians]

        """

        # A14 - De Pury & Farquhar
        decl = -23.4 * (pi / 180.) * cos(2. * pi * (doy + 10) / 365.)

        return decl

    def calculate_eqn_of_time(self, gamma):

        """
        Correction for the difference between solar time and the clock
        time - De Pury & Farquhar (1997): eq A17, also see Spencer
        (1971)

        Arguments:
        ----------
        gamma: float
            fractional year [radians]

        Returns:
        ---------
        the equation of time [minutes]

        """

        et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 *
              cos(2. * gamma) - 9.731 * sin(gamma))

        return et

    def calculate_solar_noon(self, et):

        """
        Calculation solar noon, based on De Pury & Farquhar (1997), eq
        A16

        Arguments:
        ----------
        et: float
            equation of time [radians]

        Returns:
        ---------
        t0: float
            solar noon [hours]

        """

        # all international standard meridians are multiples of 15deg
        # east/west of greenwich
        Ls = round(self.lon / 15.) * 15.
        t0 = 12. + (4. * (Ls - self.lon) - et) / 60.

        return t0

    def calculate_hour_angle(self, t, t0):

        """
        Calculation solar noon - De Pury & Farquhar (1997): eq A15

        Arguments:
        ----------
        t: float
            hour of the day [hours]

        t0: float
            solar noon [hours]

        Returns:
        ---------
        the hour angle [radians]

        """

        return pi * (t - t0) / 12.

    def calculate_solar_geometry(self, doy):

        """
        The solar zenith angle is the angle between the zenith and the
        centre of the sun's disc. The solar elevation angle is the
        altitude of the sun, the angle between the horizon and the
        centre of the sun's disc. Since these two angles are
        complementary, the cosine of either one of them equals the sine
        of the other, i.e. cos theta = sin beta. I will use cos_zen
        throughout code for simplicity.

        Arguments:
        ----------
        doy: float
            day of the year

        Returns:
        -----------
        cos_zenith: array
            diurnal time course of cosinus zenithal [degrees]

        """

        # declare empty array
        cos_zenith = np.zeros(48)

        for i in np.arange(1, 48+1):

            hod = i / 2.  # need to convert 30 min data, 0-47 to 0-23.5

            # sun's position info
            gamma = self.day_angle(doy)
            dec = self.calculate_solar_declination(doy, gamma)
            et = self.calculate_eqn_of_time(gamma)
            t0 = self.calculate_solar_noon(et)
            h = self.calculate_hour_angle(hod, t0)
            rlat = self.lat * pi / 180.

            # A13 - De Pury & Farquhar (1997)
            sin_beta = sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(h)
            cos_zenith[i-1] = sin_beta  # cos theta = sin beta

            if cos_zenith[i-1] > 1.:
                cos_zenith[i-1] = 1.

            elif cos_zenith[i-1] < 0.:
                cos_zenith[i-1] = 0.

        return cos_zenith

    def calc_extra_terrestrial_rad(self, doy, cos_zen):

        """
        Solar radiation incident outside the earth's atmosphere, e.g.
        extra-terrestrial radiation. The value varies a little with the
        earth's orbit. Using formula (eq 1) from Spitters et al. (1986).

        Arguments:
        ----------
        doy: float
            day of year

        cos_zen: float
            cosine of zenith angle [radians]

        Returns:
        --------
        So: float
            solar radiation normal to the sun's bean outside the Earth's
            atmosphere [J m-2 s-1]

        """

        # Solar constant (J m-2 s-1)
        Sc = 1362.

        if cos_zen > 0.:
            # trig funcs are cofuncs of each other
            # sin(x) = cos(90-x) and cos(x) = sin(90-x).
            S0 = Sc * (1. + 0.033 * cos(float(doy) / 365. * 2. * pi)) * cos_zen

        else:
            S0 = 0.

        return S0

    def spitters(self, doy, par_day, cos_zenith):

        """
        Spitters algorithm to estimate the diffuse component from the
        total daily incident radiation - Spitters et al. (1986): Eq.
        2a-d

        Arguments:
        ----------
        doy: int
            day of year

        par_day: float
            daily total photosynthetically active radiation
            [umol m-2 d-1]

        cos_zenith: float
            cosine of zenith angle [radians]

        Returns:
        -------
        diffuse: float
            diffuse component of incoming radiation [unitless]

        """

        # Calculate extra-terrestrial radiation
        S0 = 0.

        for i in range(48):

            S0 += (self.calc_extra_terrestrial_rad(doy, cos_zenith[i])
                   * self.SEC_2_HLFHR * self.J_TO_MJ)

        # atmospheric transmisivity
        tau = (par_day / self.SW_2_PAR_MJ) / S0

        # Spitter's formula (Eq. 2a-d)
        if tau < 0.07:
            diffuse_frac = 1.

        elif tau < 0.35:
            diffuse_frac = 1. - 2.3 * (tau - 0.07) ** 2.

        elif tau < 0.75:
            diffuse_frac = 1.33 - 1.46 * tau

        else:
            diffuse_frac = 0.23

        return diffuse_frac

    def estimate_diurnal_par(self, par_day, doy):

        """
        Calculate daily course of incident PAR from daily totals using a
        routine modified from MAESTRA
        (http://maespa.github.io/manual.html)

        Arguments:
        ----------
        par_day: float
            daily total photosynthetically active radiation
            [umol m-2 d-1]

        doy: int
            day of the year

        Returns:
        --------
        par: array
            diurnal course of the par [umol m-2 s-1]

        """

        # declare the empty arrays
        cos_bm = np.zeros(48)
        cos_df = np.zeros(48)
        par = np.zeros(48)

        # Transmissivity of the atmosphere
        tau = 0.76

        # available light?
        cos_zenith = self.calculate_solar_geometry(doy)
        diffuse_frac = self.spitters(doy, par_day, cos_zenith)
        direct_frac = 1. - diffuse_frac

        # daily total beam PAR (umol m-2 d-1)
        beam_rad = par_day * direct_frac

        # daily total diffuse PAR (umol m-2 d-1)
        diffuse_rad = par_day * diffuse_frac

        # initialisation
        sum_bm = 0.
        sum_df = 0.

        for i in np.arange(48):

            if cos_zenith[i] > 0.:
                zenith = acos(cos_zenith[i])

                # set beam = 0. for zenith > 80 degrees
                if zenith < 80. * pi / 180.:
                    cos_bm[i] = cos_zenith[i] * tau ** (1. / cos_zenith[i])

                else:
                    cos_bm[i] = 0.

                cos_df[i] = cos_zenith[i]
                sum_bm += cos_bm[i]
                sum_df += cos_df[i]

        for i in np.arange(48):

            if sum_bm > 0.:
                rdbm = beam_rad * cos_bm[i] / sum_bm

            else:
                rdbm = 0.

            if sum_df > 0.:
                rddf = diffuse_rad * cos_df[i] / sum_df

            else:
                rddf = 0.

            par[i] = (rddf + rdbm)  # umol m-2 s-1

        return par

    def estimate_diurnal_vpd(self, vpd15_prev, vpd09, vpd15, vpd09_next):

        """
        Interpolate VPD between 9am and 3pm values to generate diurnal
        VPD following the method of Haverd et al. (2013). This seems
        reasonable, vapour pressure plotted against time of day often
        does not reveal consistent patterns, with small fluctuations
        (Kimball & Bellamy (1986)).

        Arguments:
        ----------
        vpd15_prev: float
            vpd at 3 pm the previous day [kPa]

        vpd09: float
            expected vpd at 9 am [kPa]

        vpd15: float
            expected vpd at 3 pm [kPa]

        vpd09_next: float
            expected vpd at 9 am the next day [kPa]


        Returns:
        --------
        vpd: array
            diurnal course of the vapor pressure deficit [kPa]

        """

        # declare empty array
        vpd = np.zeros(48)

        # number of hours gap, i.e. 3pm to 9am the next day
        gap = 18.

        for i in np.arange(1, 48+1):

            hour = i / 2.  # convert 30 min data, i.e. 0-47 to 0-23.5

            if hour <= 9.:
                vpd[i-1] = (vpd15_prev + (vpd09 - vpd15_prev) * (9. + hour) /
                            gap)

            elif hour > 9. and hour <= 15.:
                vpd[i-1] = vpd09 + (vpd15 - vpd09) * (hour - 9.) / (15. - 9.)

            elif hour > 15.:
                vpd[i-1] = vpd15 + (vpd09_next - vpd15) * (hour - 15.) / gap

        return vpd

    def disaggregate_rainfall(self, rain_day):

        """
        Assign daily precip total to hours of the day, following MAESTRA
        (http://maespa.github.io/manual.html), which follows algorithm
        from GRAECO (cf. Loustau et al. (2001)).

        Arguments:
        ----------
        rain_day: float
            daily rainfall total [mm]

        Returns:
        ----------
        rain: array
            diurnal time course of rainfall [mm]

        """

        # declare empty array
        rain = np.zeros(48)

        # All rain falls in a couple hours for light storms (<1.5 mm)
        if rain_day <= 1.5:
            hour_index = np.random.randint(low=0, high=47)
            rain[hour_index] = rain_day

        # All rain falls in 24 hours for storms >46 mm
        elif rain_day > 46.:

            for i in np.arange(48):

                rain[i] = rain_day / 48.

        # Aim if for all rain to fall at ~1mm/hour at a random time of
        # the day. If we generate the same random number, then we
        # increase rainfall for this hour
        else:
            num_hrs_with_rain = int(rain_day)
            rate = rain_day / float(num_hrs_with_rain)

            # sample without replacement
            random_hours = np.random.randint(low=0, high=47,
                                             size=num_hrs_with_rain)

            for i in np.arange(num_hrs_with_rain):

                rain[random_hours[i]] += rate

        return rain

    def calc_day_length(self, doy, yr_days):

        """
        Daylength in hours - Leuning et al. (1995): A4, A5 and A6,

        Arguments:
        ----------
        doy: int
            day of year, 1=jan 1

        yr_days: int
            number of days in a year, 365 or 366

        Returns:
        --------
        daylength [hrs]

        """

        deg2rad = pi / 180.
        rlat = self.lat * deg2rad
        sindec = -sin(23.5 * deg2rad) * cos(2. * pi * (doy + 10.) / yr_days)
        a = sin(rlat) * sindec
        b = cos(rlat) * cos(asin(sindec))

        return 12. * (1. + (2. / pi) * asin(a / b))

    def estimate_diurnal_temp(self, doy, tmin, tmax):

        """
        Calculate diurnal temperature following Parton & Logan (1981).
        The day is divided into two segments and using a truncated sine
        wave in the daylight and an exponential decrease in temperature
        at night.

        Arguments:
        ----------
        doy: int
            day of the year

        tmin: float
            day minimum temp [deg C]

        tmax: float
            day maximum temp [deg C]

        Returns:
        ----------
        tday: array
            diurnal time course of the temperature [deg C]

        """

        # declare empty array
        tday = np.zeros(48)

        # 1.5 m air temperature from Parton & Logan (1981), table 1
        a = 1.86
        b = 2.2  # nighttime coeffcient
        c = -0.17  # lag of the min temp from the time of runrise

        day_length = self.calc_day_length(doy, 365)
        night_length = 24 - day_length
        sunrise = 12. - day_length / 2. + c
        sunset = 12. + day_length / 2.

        # temperature at sunset
        m = sunset - sunrise + c
        tset = (tmax - tmin) * sin(pi * m / (day_length + 2. * a)) + tmin

        for i in np.arange(1, 48+1):

            hour = i / 2.  # convert 30 min data, 0-47 to 0-23.5

            # hour - time of the min temperature (accounting for lag)
            m = hour - sunrise + c

            if hour >= sunrise and hour <= sunset:
                tday[i-1] = (tmin + (tmax - tmin) * sin((pi * m) /
                             (day_length + 2. * a)))

            elif hour > sunset:
                n = hour - sunset

            elif hour < sunrise:
                n = (24. + hour) - sunset
                d = (tset - tmin) / (exp(b) - 1.)

                # includes missing displacement to allow T to reach Tmin,
                # this  removes a discontinuity in the original Parton &
                # Logan eq. See Kimball & Bellamy (1986)
                tday[i-1] = ((tmin - d) + (tset - tmin - d) *
                             exp(-b * n / (night_length + c)))

        return tday
