#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Get CRU climate data for a lat, lon point (e.g. vpd, tair, AI, MAP). The
lat, lon coordinates can automatically be retrieved for a given FLUXNET
site. The climatology period is hardwired to 1972 - 2001, but that can
easily be tweaked.

The following CRU data are needed:
cru_ts4.03.1901.2018.vap.dat.nc, cru_ts4.03.1901.2018.tmp.dat.nc,
cru_ts4.03.1901.2018.tmn.dat.nc, cru_ts4.03.1901.2018.tmx.dat.nc,
cru_ts4.03.1901.2018.pre.dat.nc, cru_ts4.03.1901.2018.pet.dat.nc

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Harris, I. P. D. J., et al. "Updated high‐resolution grids of monthly
  climatic observations–the cru_ts3. 10 Dataset." International journal
  of climatology 34.3 (2014): 623-642.
* Monteith, J. L., & Unsworth, M. H. (1990). Principles of environmental
  physics. Arnold. SE, London, UK.

"""


__title__ = "CRU climate info of a specific grid point (lat, lon)"
__author__ = "Manon E. B. Sabot"
__version__ = "2.0 (20.05.2019)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

import warnings  # ignore these warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# import general modules
import argparse  # to read user input
import os  # check for files, paths
import sys  # check for version on the system
import calendar  # time info
from calendar import monthrange  # time info
import datetime  # time info
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files
import xarray as xr  # read netcdf

# own modules
try:
    from general_utils import read_csv, site_info  # read files, geoloc
    from calculate_solar_geometry import cos_zenith, top_atm_rad  # geom
    from calculate_solar_geometry import cloud_correction, composite_LAI

except (ImportError, ModuleNotFoundError):
    from TractLSM.Utils.general_utils import read_csv, site_info
    from TractLSM.Utils.calculate_solar_geometry import cos_zenith, top_atm_rad
    from TractLSM.Utils.calculate_solar_geometry import cloud_correction
    from TractLSM.Utils.calculate_solar_geometry import composite_LAI


# ======================================================================

def main(site, VPD=True, Tair=True, Txx=True, LAI=False, main_dir=None):

    """
    Main function: inserts climate info (on average or on extreme
                   climate) in the site info file and returns the
                   associated Tair and VPD.

    Arguments:
    ----------
    site: string
        name of the site which's CRU climate should be retrieved

    VPD: boolean
        if True

    Tair: boolean
        if True

    Txx: boolean
        if True

    LAI: boolean
        if True

    main_dir: string
        name

    Returns:
    --------
    var1: float
        vpd associated with the average climatological air temperature
        of the growing season or with the 90th percentile of the maximum
        monthly temperatures [kPa]

    var2: float
        average climatological air temperature of the growing season or
        90th percentile of the maximum monthly temperatures [deg C]

    Those, as well as the growing season, composite LAI, aridity index
    (AI), and mean annual precipitation (MAP [mm y-1]) are also added to
    the site info file.

    """

    if Tair and Txx:
        Tair = False

    if main_dir is None:
        main_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

    while 'src' in main_dir:
        main_dir = os.path.dirname(main_dir)

    fluxsites_folder = os.path.join(os.path.join(main_dir, 'input'),
                                    'fluxsites')

    write_csv = False
    vpd = pd.np.nan
    tair = pd.np.nan
    txx = pd.np.nan

    if any([VPD, Tair, Txx]):
        info = read_csv(os.path.join(fluxsites_folder, 'info_sites.csv'),
                        drop_units=False)

        if ('CRU cVPDxx', '[kPa]') not in info.columns.tolist():
            write_csv = True
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('CRU cVPDxx', '[kPa]'), value=pd.np.nan)

        if ('CRU cTxx', '[deg C]') not in info.columns.tolist():
            write_csv = True
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('CRU cTxx', '[deg C]'), value=pd.np.nan)

        if ('CRU cVPD', '[kPa]') not in info.columns.tolist():
            write_csv = True
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('CRU cVPD', '[kPa]'), value=pd.np.nan)

        if ('CRU cTair', '[deg C]') not in info.columns.tolist():
            write_csv = True
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('CRU cTair', '[deg C]'), value=pd.np.nan)

        if ('CRU AI', '[-]') not in info.columns.tolist():
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('CRU AI', '[-]'), value=pd.np.nan)

        if ('CRU MAP', '[mm y-1]') not in info.columns.tolist():
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('CRU MAP', '[mm y-1]'), value=pd.np.nan)

        if ('Composite LAI', '[m2 m-2]') not in info.columns.tolist():
            write_csv = True
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('Composite LAI', '[m2 m-2]'), value=pd.np.nan)

        if ('CRU growing season', '[-]') not in info.columns.tolist():
            write_csv = True
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+1,
                        column=('CRU growing season', '[-]'), value=pd.np.nan)

        columns = info.columns.values
        info.columns = info.columns.droplevel(level=1)
        info.index = info['Site']

        if site is not None:
            sites = [site]

        else:
            sites = info.index.tolist()
            sites = [str(e) for e in sites if str(e) != str(pd.np.NaN)]

        for s in sites:

            site_code = info.loc[s, 'Fluxnet Code']
            lat, lon = site_info(site_code, s)

            if VPD:  # retrieve site mean growing season VPD

                if Tair:
                    vpd = info.loc[s, 'CRU cVPD']
                    tair = info.loc[s, 'CRU cTair']

                    if ((pd.isna(vpd) or pd.isna(tair)) or (pd.isna(vpd) and
                                                            pd.isna(tair))):
                        write_csv = True

                        if pd.isna(info.loc[s, 'CRU growing season']):
                            start, stop = growing_season(fluxsites_folder, lat,
                                                         lon)

                        else:
                            start, stop = \
                                info.loc[s, 'CRU growing season'].split('-')
                            start = list(calendar.month_abbr).index(start)
                            stop = list(calendar.month_abbr).index(stop)

                        # months_2_keep = np.arange(start, stop + 1)
                        months_2_keep = np.arange(4, 11)  # Apr-Nov set
                        vpd, tair = growing_season_mean(fluxsites_folder, lat,
                                                        lon, months_2_keep)

                if Txx:
                    vpd = info.loc[s, 'CRU cVPDxx']
                    txx = info.loc[s, 'CRU cTxx']

                    if ((pd.isna(vpd) or pd.isna(txx)) or (pd.isna(vpd) and
                                                           pd.isna(txx))):
                        write_csv = True
                        vpd, txx = extreme(fluxsites_folder, lat, lon)

                else:
                    vpd = info.loc[s, 'CRU cVPD']

                    if pd.isna(vpd):
                        write_csv = True

                        if pd.isna(info.loc[s, 'CRU growing season']):
                            start, stop = growing_season(fluxsites_folder, lat,
                                                         lon)

                        else:
                            start, stop = \
                                info.loc[s, 'CRU growing season'].split('-')
                            start = list(calendar.month_abbr).index(start)
                            stop = list(calendar.month_abbr).index(stop)

                        # months_2_keep = np.arange(start, stop + 1)
                        months_2_keep = np.arange(4, 11)  # Apr-Nov set
                        vpd, __ = growing_season_mean(fluxsites_folder, lat,
                                                      lon, months_2_keep)

            else:  # retrieve site mean growing season Tair, but not VPD
                if Tair:
                    tair = info.loc[s, 'CRU cTair']

                    if pd.isna(tair):
                        write_csv = True

                        if pd.isna(info.loc[s, 'CRU growing season']):
                            start, stop = growing_season(fluxsites_folder, lat,
                                                         lon)

                        else:
                            start, stop = \
                                info.loc[s, 'CRU growing season'].split('-')
                            start = list(calendar.month_abbr).index(start)
                            stop = list(calendar.month_abbr).index(stop)

                        # months_2_keep = np.arange(start, stop + 1)
                        months_2_keep = np.arange(4, 11)  # Apr-Nov set
                        tair = growing_season_Tair(fluxsites_folder, lat, lon,
                                                   months_2_keep)

                if Txx:
                    txx = info.loc[s, 'CRU cTxx']

                    if pd.isna(txx):
                        write_csv = True
                        txx = Txx(fluxsites_folder, lat, lon)

            if write_csv and VPD:
                if not Txx:
                    info.loc[s, 'CRU cVPD'] = vpd

                else:
                    info.loc[s, 'CRU cVPDxx'] = vpd

            if write_csv and Tair:
                info.loc[s, 'CRU cTair'] = tair

            if write_csv and Txx:
                info.loc[s, 'CRU cTxx'] = txx

            # aridity index and MAP
            if (pd.isna(info.loc[s, 'CRU AI']) or
               pd.isna(info.loc[s, 'CRU MAP'])):
                write_csv = True
                info.loc[s, 'CRU AI'], info.loc[s, 'CRU MAP'] = \
                    aridity(fluxsites_folder, lat, lon)

            # growing season
            if pd.isna(info.loc[s, 'CRU growing season']):
                write_csv = True

                try:
                    info.loc[s, 'CRU growing season'] = '%s-%s' % \
                        (calendar.month_abbr[start], calendar.month_abbr[stop])

                except UnboundLocalError:
                    start, stop = growing_season(fluxsites_folder, lat, lon)
                    info.loc[s, 'CRU growing season'] = '%s-%s' % \
                        (calendar.month_abbr[start], calendar.month_abbr[stop])

            # reset vpd, tair, txx for other sites to be written as well
            if site is None and s != sites[-1]:
                write_csv = False
                vpd = pd.np.nan
                tair = pd.np.nan
                txx = pd.np.nan

    if write_csv:
        info.columns = pd.MultiIndex.from_tuples(columns)
        info.to_csv(os.path.join(fluxsites_folder, 'info_sites.csv'),
                    index=False, na_rep='', encoding='utf-8')

    write_csv = False
    lai = pd.np.nan

    if LAI:
        info = read_csv(os.path.join(fluxsites_folder, 'info_sites.csv'),
                        drop_units=False)

        if ('Composite LAI', '[m2 m-2]') not in info.columns.tolist():
            write_csv = True
            info.insert(loc=info.columns.get_loc(('Koeppen Climate', '[-]'))+2,
                        column=('Composite LAI', '[m2 m-2]'), value=pd.np.nan)

        columns = info.columns.values
        info.columns = info.columns.droplevel(level=1)
        info.index = info['Site']

        if site is not None:
            sites = [site]

        else:
            sites = info.index.tolist()
            sites = [str(e) for e in sites if str(e) != str(pd.np.NaN)]

        for s in sites:

            site_code = info.loc[s, 'Fluxnet Code']
            lat, lon = site_info(site_code, s)
            lai = info.loc[s, 'Composite LAI']

            if pd.isna(lai):
                write_csv = True

                if pd.isna(info.loc[s, 'CRU growing season']):
                    start, stop = growing_season(fluxsites_folder, lat, lon)

                else:
                    start, stop = info.loc[s, 'CRU growing season'].split('-')
                    start = list(calendar.month_abbr).index(start)
                    stop = list(calendar.month_abbr).index(stop)
                    info.loc[s, 'CRU growing season'] = '%s-%s' % \
                        (calendar.month_abbr[start], calendar.month_abbr[stop])

                # months_2_keep = np.arange(start, stop + 1)
                months_2_keep = np.arange(4, 11)  # Apr-Nov set
                lai = growing_season_LAI(info.loc[s], fluxsites_folder, lat,
                                         lon, months_2_keep)
                info.loc[s, 'Composite LAI'] = lai

            # reset vpd, tair, txx for other sites to be written as well
            if site is None and s != sites[-1]:
                write_csv = False
                lai = pd.np.nan

    if write_csv:
        info.columns = pd.MultiIndex.from_tuples(columns)
        info.to_csv(os.path.join(fluxsites_folder, 'info_sites.csv'),
                    index=False, na_rep='', encoding='utf-8')

    if site is not None:
        if VPD and Tair:

            return vpd, tair

        elif VPD and Txx:

            return vpd, txx

        elif VPD:

            return vpd, pd.np.nan

        elif Tair:

            return pd.np.nan, tair

        elif Txx:

            return pd.np.nan, txx

        elif LAI:

            return lai


# ======================================================================

# ~~~ Other functions are defined here ~~~

def CRU(fluxsites_folder, lat, lon, var='Tair'):

    """
    Accesses the relevant CRU dataset and converts units appropriately.

    Arguments:
    ----------
    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    var: string
        which CRU dataset to open?

    Returns:
    --------
    Dataframe containing the site-level data.

    """

    # read netcdf file
    CRU_folder = os.path.join(fluxsites_folder, 'CRU')

    if var == 'Tair':
        filename = 'cru_ts4.03.1901.2018.tmp.dat.nc'
        column = 'Tair'

    if var == 'Tmin':
        filename = 'cru_ts4.03.1901.2018.tmn.dat.nc'
        column = 'Tmin'

    if var == 'Tmax':
        filename = 'cru_ts4.03.1901.2018.tmx.dat.nc'
        column = 'Tmax'

    if var == 'VPD':
        filename = 'cru_ts4.03.1901.2018.vap.dat.nc'
        column = 'Vap'

    if var == 'PET':
        filename = 'cru_ts4.03.1901.2018.pet.dat.nc'
        column = 'pet'

    if var == 'PTT':
        filename = 'cru_ts4.03.1901.2018.pre.dat.nc'
        column = 'pre'

    if var == 'CLD':
        filename = 'cru_ts4.03.1901.2018.cld.dat.nc'
        column = 'cld'

    ds = xr.open_dataset(os.path.join(CRU_folder, filename))

    # retrieve site var
    lats = np.unique(ds['lat'])
    lons = np.unique(ds['lon'])
    lat = lats[np.argmin(abs(lats - lat))]
    lon = lons[np.argmin(abs(lons - lon))]
    df = ds.sel(lat=lat, lon=lon).to_dataframe()
    df.drop(columns=['lat', 'lon'], inplace=True)

    try:
        df.drop(columns=['stn'], inplace=True)

    except KeyError:
        pass

    df.columns = [column]

    if var == 'VPD':
        df['Vap'] *= 0.1  # hPa to kPa

    if var == 'PET':
        Ndpm = [monthrange(i, j)[1] for i in
                range(df.index.year[0], df.index.year[-1]+1) for j in
                range(1, 13)]  # days per month
        df['pet'] *= Ndpm  # mm/day to mm/month

    if var == 'CLD':
        df['cld'] /= 100.  # now between 0-1

    return df


def growing_season(fluxsites_folder, lat, lon):

    """
    Restricts the data to the "growing season", defined as the months
    between which the minimum temperature never drops below zero.

    Arguments:
    ----------
    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    Returns:
    --------
    The first and last month indexes for the 'growing season".

    """

    df = CRU(fluxsites_folder, lat, lon, var='Tmin')  # Tmin data

    # 1972 - 2001 period
    df = df.loc[np.logical_and(df.index.year > 1971, df.index.year < 2002)]

    to_drop = np.unique(df.index[df['Tmin'] <= 0.].month.tolist())
    remain = [e for e in np.unique(df.index.month.tolist())
              if e not in to_drop]

    # remaining months must be consecutive (e.g. no Tmin < 0 in June)
    diff_remain = [i - j for j, i in zip(remain, remain[1:])]

    # deal with non-consecutive months
    if not ((len(np.unique(diff_remain)) == 1) and (diff_remain[0] == 1)):
        idx = [i + 1 for i, e in enumerate(diff_remain) if e != 1]
        split = [remain[i:j] for i, j in zip([0] + idx,
                 idx + ([len(remain)] if idx[-1] != len(remain) else []))]
        remain = max(split, key=len)

    return [remain[0], remain[-1]]


def growing_season_mean(fluxsites_folder, lat, lon, growing_season):

    """
    Calculates the average air temperature and vapour pressure deficit
    associated with the growing season.

    Arguments:
    ----------
    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    growing_season: array
        index list for the 'growing season"

    Returns:
    --------
    vpd: float
        average vapour pressure deficit associated to tair for the
        growing season [kPa]

    tair: float
        average air temperature for the growing season [deg C]

    """

    df1 = CRU(fluxsites_folder, lat, lon)  # Tair data
    df2 = CRU(fluxsites_folder, lat, lon, var='VPD')  # VPD data

    # 1972 - 2001 period
    df1 = df1.loc[np.logical_and(df1.index.year > 1971, df1.index.year < 2002)]
    df2 = df2.loc[np.logical_and(df2.index.year > 1971, df2.index.year < 2002)]

    # sat vapor pressure (Tetens eq., as in Monteith & Unsworth, 1990)
    es = 0.61078 * np.exp(17.27 * df1['Tair'] / (df1['Tair'] + 237.3))

    # update vpd to go from vap to vpd
    df2['Vap'] = es - df2['Vap']
    df2.columns = ['VPD']

    # average the monthly values across years & select growing season
    Tair_mean = (df1.groupby(df1.index.month).mean()).loc[growing_season]
    VPD_mean = (df2.groupby(df2.index.month).mean()).loc[growing_season]

    # get yearly average of multi-year monthly values
    tair = Tair_mean.mean().values[0]
    vpd = VPD_mean.mean().values[0]

    return vpd, tair


def extreme(fluxsites_folder, lat, lon):

    """
    Calculates the 90th percentile of maximum monthly air temperatures
    and associated vapour pressure deficit.

    Arguments:
    ----------
    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    Returns:
    --------
    vpd: float
        vapour pressure deficit associated to the 90th percentile of
        maximum monthly tair [kPa]

    txx: float
        90th percentile of maximum monthly air temperature [deg C]

    """

    df1 = CRU(fluxsites_folder, lat, lon, var='VPD')  # VPD data
    df2 = CRU(fluxsites_folder, lat, lon, var='Tmax')  # Tmax data

    # 1972 - 2001 period
    df1 = df1.loc[np.logical_and(df1.index.year > 1971, df1.index.year < 2002)]
    df2 = df2.loc[np.logical_and(df2.index.year > 1971, df2.index.year < 2002)]

    # sat vapor pressure (Tetens eq., as in Monteith & Unsworth, 1990)
    es = 0.61078 * np.exp(17.27 * df2['Tmax'] / (df2['Tmax'] + 237.3))

    # update vpd to go from vap to vpd
    df1['Vap'] = es - df1['Vap']
    df1.columns = ['VPD']

    # select 90th percentile of max temperatures & associated VPD
    N = int(len(df2) / 10)
    df1 = df1[(df2 >= df2.nlargest(N, 'Tmax').min())['Tmax']]  # VPD
    df2 = df2.nlargest(N, 'Tmax')  # N largest T values
    df1 = df1.loc[df2.nsmallest(1, 'Tmax').index]  # VPD
    df2 = df2.nsmallest(1, 'Tmax')  # actual 90th percentile value
    vpd = df1.values[0][0]
    txx = df2.values[0][0]

    return vpd, txx


def growing_season_Tair(fluxsites_folder, lat, lon, growing_season):

    """
    Calculates the average air temperature associated with the growing
    season.

    Arguments:
    ----------
    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    growing_season: array
        index list for the 'growing season"

    Returns:
    --------
    tair: float
        average air temperature for the growing season [deg C]

    """

    df = CRU(fluxsites_folder, lat, lon)  # Tair data

    # 1972 - 2001 period
    df = df.loc[np.logical_and(df.index.year > 1971, df.index.year < 2002)]

    # average the monthly values across years & select growing season
    Tair_mean = (df.groupby(df.index.month).mean()).loc[growing_season]

    return Tair_mean.mean().values[0]


def Txx(fluxsites_folder, lat, lon):

    """
    Calculates the 90th percentile of maximum monthly air temperatures.

    Arguments:
    ----------
    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    Returns:
    --------
    txx: float
        90th percentile of maximum monthly air temperature [deg C]

    """

    df = CRU(fluxsites_folder, lat, lon, var='Tmax')  # Tmax data

    # 1972 - 2001 period
    df = df.loc[np.logical_and(df.index.year > 1971, df.index.year < 2002)]

    # select 90th percentile of max temperatures
    N = int(len(df) / 10)
    df = df.nlargest(N, 'Tmax')  # N largest T values
    df = df.nsmallest(1, 'Tmax')  # actual 90th percentile value

    return df.values[0][0]


def aridity(fluxsites_folder, lat, lon):

    """
    Calculates the aridity index and mean annual precipitation.

    Arguments:
    ----------
    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    Returns:
    --------
    AI: float
        aridity index [unitless]

    MAP: float
        mean annual precipitation [mm y-1]

    """

    df1 = CRU(fluxsites_folder, lat, lon, var='PET')  # PET data
    df2 = CRU(fluxsites_folder, lat, lon, var='PTT')  # precip data

    # 1972 - 2001 period
    df1 = df1.loc[np.logical_and(df1.index.year > 1971, df1.index.year < 2002)]
    df2 = df2.loc[np.logical_and(df2.index.year > 1971, df2.index.year < 2002)]

    # AI
    yPET = df1.groupby(df1.index.year).sum()
    yPTT = df2.groupby(df2.index.year).sum()
    AI = yPTT.values / yPET.values

    return AI.mean(), yPTT.values.mean()


def growing_season_LAI(ds, fluxsites_folder, lat, lon, growing_season):

    """
    Calculates the average LAI of the weighted sun/shade LAI over the
    growing season.

    Arguments:
    ----------
    ds: pandas series
        site info for which to calculate the LAI

    fluxsites_folder: string
        path to the folder where the site info should be stored

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    growing_season: array
        index list for the 'growing season"

    Returns:
    --------
    LAI: float
        idealised composite LAI [m2 m-2]

    """

    # LAI climatologies
    site = ds['Site']
    df1, __ = read_csv(os.path.join(os.path.join(fluxsites_folder, 'LAI'),
                                    '%s_lai_climatology.csv' % (site)))

    # CRU climate data
    df2 = CRU(fluxsites_folder, lat, lon, var='CLD')  # cloud cover data

    # 1972 - 2001 period
    df2 = df2.loc[np.logical_and(df2.index.year > 1971, df2.index.year < 2002)]

    # April - November, hardcoded growing season
    df2 = df2.loc[np.logical_and(df2.index.month > 3, df2.index.month < 11)]

    # add year, month column
    df2['year'] = df2.index.year
    df2['month'] = df2.index.month

    # growing season doys, hardcoded between Apr-Nov for a 365 day year
    start = 91
    end = 305
    years = np.repeat(np.unique(df2.index.year),
                      len(np.repeat(np.arange(start, end), 24)))
    doy = np.tile(np.repeat(np.arange(start, end), 24),
                  len(np.unique(df2.index.year)))
    months = [(datetime.datetime(years[i], 1, 1) +
               datetime.timedelta(doy[i] - 1)).month
              for i in range(len(years))]
    hod = np.tile(np.tile(np.arange(0, 24), end - start),
                  len(np.unique(df2.index.year)))

    # cosinus zenithal
    coszen = cos_zenith(doy, hod, lat, lon)

    # extra-terrestrial incident rad, eq 1 Spitters et al. (1986)
    Sxt = top_atm_rad(doy, coszen)  # umol m-2 s-1

    # put Sxt and coszen into a dataframe
    df3 = pd.DataFrame(np.asarray([years, months, doy, hod, coszen, Sxt]).T,
                       columns=['year', 'month', 'doy', 'hod', 'coszen',
                                'Sxt'])

    # remove leap year day in March
    df3 = df3[df3['month'] != 3.]

    # append cloud cover, temperature, VPD, LAI to df5
    df2 = df2.reset_index().drop(columns=['time'])
    df3 = df3.merge(df2, how='inner', on=['year', 'month'])
    df3 = df3.merge(df1, how='inner', on=['doy'])

    # correct the radiation using the cloud cover
    df3['PPFD'] = df3['Sxt']
    df3['PPFD'][df3['Sxt'] > 0.] = \
        cloud_correction(df3['Sxt'][df3['Sxt'] > 0.],
                         df3['cld'][df3['Sxt'] > 0.], lat)

    # declare the new weighted composite LAI variable
    df3['wcLAI'] = df3['LAI']
    df3['wcLAI'][df3['PPFD'] <= 0.] = 0.

    # restrict the timeseries (for speed)
    doys = np.asarray(df3['doy'][df3['PPFD'] > 0.].values)
    cos_zen = np.asarray(df3['coszen'][df3['PPFD'] > 0.].values)
    Sxts = np.asarray(df3['Sxt'][df3['PPFD'] > 0.].values)
    PPFDs = np.asarray(df3['PPFD'][df3['PPFD'] > 0.].values)
    LAIs = np.asarray(df3['LAI'][df3['PPFD'] > 0.].values)

    # scatter the LAI over the whole period
    df3['wcLAI'][df3['PPFD'] > 0.] = composite_LAI(doys, cos_zen, Sxts, PPFDs,
                                                   LAIs)
    LAI_comp = np.asarray(df3['wcLAI'].values)

    # average new composite LAI per day
    LAI_comp = LAI_comp[:(len(LAI_comp) // 24) * 24].reshape(-1, 24)
    LAI_comp = np.ma.masked_where(LAI_comp <= 0., LAI_comp)

    return np.mean(np.ma.mean(LAI_comp, axis=1))


if __name__ == "__main__":

    # define the argparse settings to read run set up file
    description = "Retrieve mean VPD and/or temperature of growing season at \
                   a specific site"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-s', '--site', type=str, help='site name (kPa)',
                        default=None)
    parser.add_argument('-VPD', '--VPD', action='store_true',
                        help='get long term VPD (kPa)')
    parser.add_argument('-Tair', '--Tair', action='store_true',
                        help='get long term Tair (C)')
    parser.add_argument('-Txx', '--Txx', action='store_true',
                        help='get long term 90th percentile Tair (C)')
    parser.add_argument('-LAI', '--LAI', action='store_true',
                        help='get site composite LAI (m2 m-2)')
    args = parser.parse_args()

    # main parent directory
    script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
    main_dir = os.path.dirname(os.path.dirname(script_dir))

    if args.site is None:
        main(args.site, VPD=args.VPD, Tair=args.Tair, Txx=args.Txx,
             LAI=args.LAI, main_dir=main_dir)

    else:
        if args.VPD:
            vpd, tair = main(args.site, VPD=args.VPD, Tair=args.Tair,
                             Txx=args.Txx, LAI=args.LAI, main_dir=main_dir)
            print(vpd, tair)

        elif args.Tair:
            vpd, tair = main(args.site, VPD=args.VPD, Tair=args.Tair,
                             Txx=args.Txx, LAI=args.LAI, main_dir=main_dir)
            print(vpd, tair)

        elif args.LAI:
            lai = main(args.site, VPD=args.VPD, Tair=args.Tair, Txx=args.Txx,
                       LAI=args.LAI, main_dir=main_dir)
            print(lai)
