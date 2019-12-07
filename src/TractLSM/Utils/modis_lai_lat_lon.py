#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Get MODIS LAI data for a lat, lon point or area around a point. The lat,
lon coordinates can automatically be retrieved for a given FLUXNET site.

Useful definitions:

LAI QC:
    https://lpdaac.usgs.gov/sites/default/files/public/
            product_documentation/mod15_user_guide.pdf
LAI SD:
    Once all 30 models are created they are applied to the MODIS data to
    yield 30 independent results. The 30 independent results are
    averaged together to yield one result for any given pixel. Standard
    deviation from the 30 results is retained in a QA layer for the end
    user to understand the amount of agreement between the independent
    models.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "LAI around specific grid point (lat, lon coordinates)"
__author__ = "Manon E. B. Sabot"
__version__ = "2.0 (16.10.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# import general modules
import os  # check for files, paths
import sys  # check for version on the system
import multiprocessing  # multiprocessing
from joblib import Parallel, delayed  # multiprocessing
import datetime  # to set the end year of retrieval to current year
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files
import requests  # request data
import json  # deal with data in json database format

# own modules
try:
    from general_utils import site_info  # site geolocation

except (ImportError, ModuleNotFoundError):
    from TractLSM.Utils.general_utils import site_info


# ======================================================================

def main(LAI_folder, fluxsites_folder, site=None, start_year=None,
         end_year=None, above_below=None, left_right=None):

    """
    Main function: reads in site-level info from the site info file,
                   then retrieves the LAI data at site-level following a
                   set of user defined specifications. If no site is
                   defined, the LAI for each of the sites present in the
                   info file will be retrieved.

    Arguments:
    ----------
    LAI_folder: string
        path to the folder where the LAI should be stored

    fluxsites_folder: string
        path to the folder where the site-level met files are stored

    site: string
        name of the site which's LAI should be retrieved

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
    A csv file with a name of the form site_filtered_scaled_Lai_500m.csv
    in the LAI_folder.

    """

    # site-level FLUXNET info
    info = pd.read_csv(os.path.join(fluxsites_folder, 'info_sites.csv'),
                       header=[0], skiprows=0)
    info.index = info['Site']
    info.drop(info.loc['[-]'], axis=0, inplace=True)

    print('Retrieving the MODIS LAI data from the NASA servers: this is long!')

    if site is None:

        sub_sites = info.index.tolist()  # index object to list object
        sub_sites = split(sub_sites, len(sub_sites) / Njobs)

        for sites in sub_sites:

            (Parallel(n_jobs=Njobs)(delayed(MODIS_LAI_csv)(LAI_folder,
                                    sites[e], info, start_year=start_year,
                                    end_year=end_year, above_below=above_below,
                                    left_right=left_right)
                                    for e in range(len(sites))))

    else:
        MODIS_LAI_csv(LAI_folder, site, info, start_year=start_year,
                      end_year=end_year, above_below=above_below,
                      left_right=left_right)

    return


# ======================================================================

# ~~~ Other functions are defined here ~~~

class retrieval_config(object):

    def __init__(self):

        """
        Default specifications for the site-level data retrieval.

        """

        # defautl MODIS info
        self.url = "https://modis.ornl.gov/rst/api/v1/"  # MODIS servers
        self.prod = 'MCD15A2H'  # MODIS/Terra+Aqua LAI/FPAR 8-Day L4
        self.band = 'Lai_500m'  # LAI band
        self.sd_band = "LaiStdDev_500m"  # LAI SD band
        self.qc_band = 'FparLai_QC'  # LAI QC band
        self.header = {'Accept': 'text/json'}  # MODIS database format

        # the user can change the year & above_below as they please
        self.start_year = 2000
        self.end_year = datetime.datetime.now().year
        self.above_below = 0.5  # km above/below
        self.left_right = 0.5  # km left/right

        return


def split(a, N):

    """
    Splits a list or array into N-roughly equal parts.

    Arguments:
    ----------
    a: list or array
        list/array to be split, can contain any data type

    N: int
        number of sub-lists/arrays the input must be split into

    Returns:
    --------
    A list of N-roughly equal parts.

    """

    integ = int(len(a) / N)
    remain = int(len(a) % N)

    splitted = [a[i * integ + min(i, remain):(i + 1) * integ +
                  min(i + 1, remain)] for i in range(N)]

    return splitted


def get_dates(url, header, prod, lat, lon, start_year, end_year):

    """
    Get the dates of the requested MODIS time series in user readable
    formats.

    Arguments:
    ----------
    url: string
        MODIS server access url

    header: dictionary
        which database format does the user want?

    prod: string
        version of the MODIS product requested

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    start_year: int
        start data retrieval that year

    end_year: int
        end data retrieval that year

    Returns:
    --------
    modis_dates: array
        MODIS-specific formatted dates

    dates: array
        calendar dates

    """

    request_url = "%s%s/dates?latitude=%s&longitude=%s" % (url, prod, str(lat),
                                                           str(lon))
    response = requests.get(request_url, headers=header)

    istart = 0
    start_date = json.loads(response.text)['dates'][istart]['modis_date']

    while int(start_date[1:5]) < start_year:

        istart += 1
        start_date = json.loads(response.text)['dates'][istart]['modis_date']

    iend = -1
    end_date = json.loads(response.text)['dates'][iend]['modis_date']

    while int(end_date[1:5]) > end_year:

        iend -= 1
        end_date = json.loads(response.text)['dates'][iend]['modis_date']

    dates = json.loads(response.text)['dates'][istart:iend]
    modis_dates = [str(e['modis_date']) for e in dates]
    dates = [str(e['calendar_date']) for e in dates]

    return modis_dates, dates


def request_data(url, prod, lat, lon, band, dates, dt1, dt2, above_below,
                 left_right):

    """
    Write the request for the MODIS data.

    Arguments:
    ----------
    url: string
        MODIS server access url

    prod: string
        version of the MODIS product requested

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    band: string
        LAI band requested

    dates: array
        MODIS-specific formatted dates

    dt1: int
        index of the first date of the request

    dt2: int
        index of the last date of the request

    above_below: float
        distance above and below selected grid point on a flat grid (km)

    left_right: float
        distance left and right of selected grid point on a flat grid
        (km)

    Returns:
    --------
    url of the request (string)

    """

    AB = str(int(round(above_below)))
    LR = str(int(round(left_right)))

    request_url = (("%s%s/subset?latitude=%s&longitude=%s&band=%s&startDate=%s"
                    + "&endDate=%s&kmAboveBelow=%s&kmLeftRight=%s") % (url,
                   prod, str(lat), str(lon), band, dates[dt1], dates[dt2], AB,
                   LR))

    return request_url


def get_data(data_array, url, header, prod, lat, lon, band, dates, dt1, dt2,
             above_below, left_right, scale=True, long_name=False):

    """
    Actually retrieve the data and append it to an existing data array

    Arguments:
    ----------
    data_array: array
        the data is appended to that array, so it could already contain
        data

    url: string
        MODIS server access url

    header: dictionary
        which database format does the user want?

    prod: string
        version of the MODIS product requested

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    band: string
        LAI band requested

    dates: array
        MODIS-specific formatted dates

    dt1: int
        index of the first date of the request

    dt2: int
        index of the last date of the request

    above_below: float
        distance above and below selected grid point on a flat grid (km)

    left_right: float
        distance left and right of selected grid point on a flat grid
        (km)

    scale: boolean
        if True, the scale info is also included

    long_name: boolean
        if True the MODIS long name corresponding to the data is also
        returned

    Returns:
    --------
    data_array: array
        MODIS LAI data

    long_name: string
        MODIS long name corresponding to the data

    """

    # Build data request url and submit request
    request_url = request_data(url, prod, lat, lon, band, dates, dt1, dt2,
                               above_below, left_right)
    response = requests.get(request_url, headers=header)

    # Loop through list of dictionaries inside the response subset key
    if scale:
        sc = float(json.loads(response.text)['scale'])

    for tstep in json.loads(response.text)['subset']:
        vals = tstep['data']

        if scale:
            vals = [i * sc for i in vals]

        data_array.append(vals)  # append data

        if long_name:
            long_name.append('%s.%s.%s.%s.%s' % (prod, tstep['modis_date'],
                             tstep['tile'], tstep['proc_date'], band))

    if long_name:
        return data_array, long_name

    return data_array


def correct_space(lai, sd, qc, above_below, left_right):

    """
    Full km footprints are always returned around a data point, no
    matter the request, with a minimum footprint of 1.5 * 1.5 km.
    Correct for this by sub-selecting the points actually requested by
    the user.

    Arguments:
    ----------
    lai: array
        lai data straight out from the MODIS servers (m2 m-2)

    sd: array
        standard deviations on the data straight out from the MODIS
        servers

    qc: array
        quality flags on the data straight out from the MODIS servers

    above_below: float
        distance above and below selected grid point on a flat grid (km)

    left_right: float
        distance left and right of selected grid point on a flat grid
        (km)

    Returns:
    --------
    lai: array
        sub-selection of lai data matching the user defined above and
        below and left and right

    sd: array
        sub-selection of standard deviations matching the user defined
        above and below and left and right

    qc: array
        sub-selection of quality flags matching the user defined above
        and below and left and right

    """

    iactual = (len(lai[0]) - 1) / 2

    if above_below < round(above_below):
        vertical = int(above_below / 0.5)

    else:
        vertical = 1

    if left_right < round(left_right):
        horizontal = int(left_right / 0.5)

    else:
        horizontal = 1

    # select the data
    if int(round(left_right)) != 0:
        ifirst = iactual - 3 * (vertical * horizontal + 1)
        ilast = iactual + 3 * (vertical * horizontal + 1)

    if int(round(left_right)) == 0:
        ifirst = iactual - 3 * (vertical) - 1
        ilast = iactual + 3 * (vertical) + 1

    ilist = []
    iter = ifirst

    while iter < ilast:
        ilist.append(np.arange(iter, iter + horizontal + vertical + 1))

        if int(round(left_right)) != 0:
            iter += horizontal + vertical + 3

        if int(round(left_right)) == 0:
            iter += horizontal + vertical + 1

    ilist = [e for sublist in ilist for e in sublist]
    lai = [[e[i] for i in ilist] for e in lai]
    sd = [[e[i] for i in ilist] for e in sd]
    qc = [[e[i] for i in ilist] for e in qc]

    return lai, sd, qc


def screen_data(lai, sd, qc, dates):

    """
    Keeps the highest quality data only.

    Arguments:
    ----------
    lai: array
        lai data (m2 m-2)

    sd: array
        standard deviations on the data

    qc: array
        quality flags on the data

    dates: array
        calendar dates for the data

    Returns:
    --------
    lai: array
        lai data (m2 m-2) screened for the highest quality

    """

    # QC screening
    qc = pd.DataFrame(qc, index=dates)

    # 'invalid' values (leading zero dropped, see Table 5 QC guide)
    mask = [0b1, 0b11, 0b100, 0b101, 0b110, 0b111, 0b1000, 0b1001, 0b1010,
            0b1011, 0b1100, 0b1101, 0b1110, 0b1111, 0b1001000, 0b1001001,
            0b1001010, 0b1001011, 0b1001100, 0b1001101, 0b1001110, 0b1001111,
            0b1101000, 0b1101001, 0b1101010, 0b1101011, 0b1101100, 0b1101101,
            0b1101110, 0b1101111, 0b10000111, 0b10000110, 0b10000101,
            0b10000100, 0b10000011, 0b10000010, 0b10000000, 0b10000001,
            0b10001111, 0b10001110, 0b10001101, 0b10001100, 0b10001011,
            0b10001010, 0b10001000, 0b10001001, 0b10010111, 0b10010110,
            0b10010101, 0b10010100, 0b10010011, 0b10010010, 0b10010000,
            0b10010001, 0b10011111, 0b10011110, 0b10011101, 0b10011100,
            0b10011011, 0b10011010, 0b10011000, 0b10011001]  # int
    mask += [255, 254, 252, 251, 250, 248]  # Fill Values
    qc.replace(mask, False)
    qc = np.asarray([True if e is not False else False for e in qc])

    # screen LAI data
    lai = pd.DataFrame(lai, index=pd.to_datetime(dates))[qc]

    # SD screening
    sd = pd.DataFrame(sd, index=[(pd.to_datetime(dates).year),
                                 (pd.to_datetime(dates).month)])
    monthly_means = lai.groupby([(pd.to_datetime(dates).year),
                                 (pd.to_datetime(dates).month)]).mean()
    iNrepeat = sd.groupby(level=[0, 1]).size()
    monthly_means = monthly_means.loc[monthly_means.index.repeat(iNrepeat)]

    # impose two arbritrary upper-values thresholds on sd
    sd.where(np.logical_and(sd < 6. * monthly_means, sd < 12.), inplace=True)
    sd.index = pd.to_datetime(dates)

    # re-screen LAI data
    lai.where(sd.notnull(), inplace=True)
    lai.where(lai > 0., inplace=True)

    return lai


def modis_screened_lai(url, header, lat, lon, prod, band, sd_band, qc_band,
                       above_below, left_right, start_year, end_year):

    """
    Retrieves the data by increments of 10 (max allowed by the MODIS
    server), sub-selects an area around a data point, and filters to
    only keep the highest quality data.

    Arguments:
    ----------
    url: string
        MODIS server access url

    header: dictionary
        which database format does the user want?

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    prod: string
        version of the MODIS product requested

    band: string
        LAI band requested

    sd_band: string
        standard deviation band requested

    qc_band: string
        quality flag band requested

    above_below: float
        distance above and below selected grid point on a flat grid (km)

    left_right: float
        distance left and right of selected grid point on a flat grid
        (km)

    start_year: int
        start data retrieval that year

    end_year: int
        end data retrieval that year

    Returns:
    --------
    dataframe containing all the subselected screened LAI data.

    """

    mdates, dates = get_dates(url, header, prod, lat, lon, start_year,
                              end_year)

    lai_data = []
    long_names = []
    sd_data = []
    qc_data = []

    # loop in increments of 10, so we don't crash the MODIS server
    for dt in range(0, len(mdates) - 10, 10):

        print(dt)

        dt1 = dt
        dt2 = dt + 9

        lai_data, long_names = get_data(lai_data, url, header, prod, lat, lon,
                                        band, mdates, dt1, dt2, above_below,
                                        left_right, long_name=long_names)
        sd_data = get_data(sd_data, url, header, prod, lat, lon, sd_band,
                           mdates, dt1, dt2, above_below, left_right)
        qc_data = get_data(qc_data, url, header, prod, lat, lon, qc_band,
                           mdates, dt1, dt2, above_below, left_right,
                           scale=False)

    # we may still have dates left over, get the last <= 10 dates
    if dates[dt2] < dates[-1]:
        lai_data, long_names = get_data(lai_data, url, header, prod, lat, lon,
                                        band, mdates, dt2 + 1, -1, above_below,
                                        left_right, scale=True,
                                        long_name=long_names)
        sd_data = get_data(sd_data, url, header, prod, lat, lon, sd_band,
                           mdates, dt2 + 1, -1, above_below, left_right)
        qc_data = get_data(qc_data, url, header, prod, lat, lon, qc_band,
                           mdates, dt2 + 1, -1, above_below, left_right,
                           scale=False)

    # if above_below and/or left_right wasn't an int, deal with it
    if (((above_below < round(above_below)) or
         (left_right < round(left_right))) or
        ((above_below < round(above_below)) and
            (left_right < round(left_right)))):
        lai_data, sd_data, qc_data = correct_space(lai_data, sd_data, qc_data,
                                                   above_below, left_right)

    if len(long_names) > len(lai_data):
        iactual = (len(lai_data[0]) - 1) / 2 + 1
        long_names = [e[iactual] for e in long_names]

    # screen data
    lai = screen_data(lai_data, sd_data, qc_data, dates)

    # average lai columns (screened) if several
    lai = lai.mean(axis=1, skipna=True)
    df = lai.reset_index(drop=True)
    df.index = mdates
    df = df.to_frame()
    df.columns = [6]

    # insert all other columns with relevant MODIS info
    df.insert(loc=0, column=5, value=[band] * len(df))

    try:
        df.insert(loc=0, column=4,
                  value=[e.split('.')[-2] for e in long_names])

    except IndexError:  # when part of the long name is missing
        df.insert(loc=0, column=4, value=long_names)

    df.insert(loc=0, column=3, value=['Lat%fLon%f' % (lat, lon)] * len(df))
    df.insert(loc=0, column=2, value=mdates)
    df.insert(loc=0, column=1, value=[prod] * len(df))
    df.insert(loc=0, column=0, value=long_names)

    return df


def MODIS_LAI_csv(LAI_folder, site, df, start_year=None, end_year=None,
                  above_below=None, left_right=None):

    """
    Builds the MODIS LAI into a dataframe and a csv file in the
    LAI_folder.

    Arguments:
    ----------
    LAI_folder: string
        path to the folder where the LAI should be stored

    site: string
        name of the site which's LAI should be retrieved

    df: dataframe
        contains site-level FLUXNET info

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
    A csv file with a name of the form site_filtered_scaled_Lai_500m.csv
    in the LAI_folder.

    """

    print(site)

    fname = '%s_filtered_scaled_Lai_500m.csv' % (os.path.join(LAI_folder,
                                                              site))

    site_code = df.loc[site, 'Fluxnet Code']
    lat, lon = site_info(site_code, site)

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
                                d.sd_band, d.qc_band, above_below, left_right,
                                start_year, end_year)
        df.to_csv(fname, header=False, index=False, na_rep='F',
                  encoding='utf-8')

    except Exception:  # any exception, still do other sites
        pass

    return


if __name__ == "__main__":

    Njobs = multiprocessing.cpu_count()

    if Njobs > 4:
        Njobs = 8

    script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
    upper_dir = os.path.dirname(os.path.dirname(script_dir))
    fluxsites_folder = os.path.join(os.path.join(upper_dir, 'input'),
                                    'fluxsites')
    LAI_folder = os.path.join(fluxsites_folder, 'LAI')
    main(LAI_folder, fluxsites_folder)
