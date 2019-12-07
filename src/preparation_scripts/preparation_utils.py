# -*- coding: utf-8 -*-

"""
General support functions

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "useful general support functions"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (18.01.2019)"
__email__ = "m.e.b.sabot@gmail.com"

#=======================================================================

# import general modules
import os  # check for files, paths
import sys  # check for files, paths
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files

try:
    import urllib2  # interrogate the FLUXNET servers for lat, lon info

except (ImportError, ModuleNotFoundError):
    import urllib

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.SPAC import Weibull_params
from TractLSM.Utils.modis_lai_lat_lon import retrieval_config
from TractLSM.Utils.modis_lai_lat_lon import modis_screened_lai
from TractLSM.Utils.met_flux_LAI_site_level import LAI_climatology


#=======================================================================

def slope(P50, which='Christoffersen'):

    """
    Finds the slope of a plant's vulnerability curve, using the water
    potential associated with a 50% decrease in hydraulic conductance as
    well as a statistical relationship, either updated from
    Christoffersen et al. (2016), which was specifically developped for
    tropical forests, or from Martin-StPaul et al. (2017), which was not
    developped with an eye to tropical forests.

    Arguments:
    ----------
    P50: float
        water potential at 50% decrease in hydraulic conductance [-MPa]

    Returns:
    --------
    sl: float
        the slope of the plant's vulnerability curve [% MPa-1]

    """

    if which == 'Christoffersen':  # (Christoffersen, pers. com.)
        sl = 65.15 * (P50) ** (-1.25)  # % MPa-1

    if which == 'Martin-StPaul':  # (Martin-StPaul et al. 2017)
        sl = 16. + np.exp(-P50) * 1092  # % MPa-1

    return sl


def Px(p, x):

    """
    Finds the leaf water potential associated with a specific x%
    decrease in hydraulic conductance, using the plant vulnerability
    curve.

    Arguments:
    ----------
    p: pandas series
        met forcings & params

    x: float
        percentage loss in hydraulic conductance

    Returns:
    --------
    Px: float
        leaf water potential [MPa] at which x% decrease in hydraulic
        conductance is observed

    """

    b, c = Weibull_params(p)  # MPa, unitless

    Px = b * ((- np.log(1 - float(x) / 100.)) ** (1. / c))  # -MPa

    return Px


def van_Genuchten_to_Clapp_Hornberger(d, m, n):

    """
    Obviously, this is to be avoided, BUT...
    Finds the b parameter from the Clapp Hornberger model, by
    equivalence between

    ((((theta - theta_hy) / (theta_sat - theta_hy)) ** (-1. / m) - 1)
     ** (1. / n))

    and

    (theta / theta_sat) ** (-b)


    i.e van Genuchten and Clapp Hornberger

    """

    bch = (np.log(((d.fc - d.pwp) / (d.theta_sat - d.pwp)) ** (-1. / m) - 1) /
           (n * (np.log(d.theta_sat) - np.log(d.fc))))

    return bch


def tetens_eq(T, RH):

    """
    Calculates vapor pressure deficit using Tetens equation, as written
    in Monteith and Unsworth (1990).

    Arguments:
    ----------
    T: float or array
        temperature [deg C]

    RH: float or array
        relative humidity [%]

    Returns:
    --------
    vpd: float or array
        vapour pressure deficit [kPa]

    """

    # saturation vapor pressure (Tetens eq.)
    es = 0.61078 * np.exp(17.27 * T / (T + 237.3))  # kPa
    ea = es * RH / 100.  # kPa

    return es - ea


def CO2_from_NOAA(df):

    """
    Access info NOAA global monthly CO2 concentrations for the time
    series present in a pandas dataset and returns it in the format
    matching the dataset's length and units.

    Arguments:
    ----------
    df: pandas dataframe
        input data

    Returns:
    --------
    CO2: array
        world monthly CO2 concentration repeated as many times as needed
        [Pa]

    """

    success = True

    years = np.unique(df['year'])

    try:
        url = "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt"

        if sys.version_info < (3, 0):  # python 2.7 or older
            response = urllib2.urlopen(url, timeout=10)

        else:
            response = urllib.request.urlopen(url, timeout=10)

    except Exception:
        success = False

        return None

    if success:
        data = response.read()

        if sys.version_info >= (3, 0):
            data = data.decode('utf-8')

        data = data.split('\n')
        data = [e.split() for e in data if '#' not in e][:-1]
        data = [float(e[3]) for e in data if any(years == int(e[0]))]

        CO2 = []

        for i in range(len(years)):

            data = [np.repeat([data[j + i * 12]],
                              len(df[np.logical_and(df['year'] == years[i],
                                                    df['month'] == j + 1)]))
                    for j in range(12)]

            data = np.asarray([item for sublist in data for item in sublist])

            CO2 = np.concatenate([CO2, data])

        return CO2 / 101.325 * 10.  # Pa


def build_LAI_climatology(fname, lat, lon, start_year=None, end_year=None,
                          above_below=None, left_right=None):

    """
    Calculates the LAI climatology based on MODIS LAI data retrieved
    from the MODIS servers if it's not already present. As MODIS
    contains NaNs (excluded) a shut-down period of 4/5 days each year
    (excluded), and the product used is with 4/8 day frequency, the
    climatology is shorter than the actual number of days in a year and
    needs to be extrapolated to 366 days.

    Arguments:
    ----------
    fname: string
        name (path) where to store the MODIS LAI data

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    start_year: int
        start data retrieval that year

    end_year: int
        end data retrieval that year

    above_below: float
        distance above and below selected grid point on a flat grid (km)

    left_right: float
        distance left and right of selected grid point on a flat grid
        (km)

    Returns:
    --------
    A pandas dataframe containing the multi-year climatology.

    """

    if not os.path.isfile(fname):

        d = retrieval_config()

        if start_year is None:
            start_year = d.start_year

        if end_year is None:
            end_year = d.end_year

        if above_below is None:
            above_below = d.above_below

        if left_right is None:
            left_right = d.left_right

        try:
            df = modis_screened_lai(d.url, d.header, lat, lon, d.prod, d.band,
                                    d.sd_band, d.qc_band, above_below,
                                    left_right, start_year, end_year)
            df.to_csv(fname, header=False, index=False, na_rep='F',
                      encoding='utf-8')

        except Exception:  # any exception, still do other sites
            print('The MODIS LAI retrieval failed')

    df = LAI_climatology(fname)

    return df
