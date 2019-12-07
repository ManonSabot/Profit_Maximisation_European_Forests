# -*- coding: utf-8 -*-

"""
General support functions

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "useful general support functions"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (18.01.2019)"
__email__ = "m.e.b.sabot@gmail.com"

# ======================================================================

# import general modules
import os  # check for paths
import sys  # check for files, version on the system
import inspect  # check for object type
import warnings  # warn the user
import importlib  # dynamic module import
import gc  # free memory when reading netcdfs
import pandas as pd  # read/write dataframes, csv files
import xarray as xr  # read netcdf

try:
    import urllib2  # interrogate the FLUXNET servers for lat, lon info

except (ImportError, ModuleNotFoundError):
    import urllib  # interrogate the FLUXNET servers for lat, lon info


# ======================================================================

def get_script_dir():

    """
    Returns a script's directory.

    """

    return os.path.dirname(os.path.realpath(sys.argv[0]))


def get_main_dir():

    """
    Returns the parent directory of a script's directory

    """

    return os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))


def retrieve_class(fname, cname=None):

    """
    Retrieves class from py file/module. If the file contains multiple
    classes, the name of the class to fetch will have to be supplied.

    Arguments:
    ----------
    fname: string
        python input filename (with path), must be in current directory
        or sub-directory(ies)

    cname: string
        name of class to retrieve

    Returns:
    --------
    cl: python class
        the one you went looking for

    """

    # transform path name into python readable module name
    if get_script_dir() in fname:
        fname = (fname.split(get_script_dir())[1].split('.py')[0]
                 .replace('/', '.').replace('\\', '.'))

    else:
        warnings.warn('The parameter file must be in the directory this ' +
                      'script is or at a lower level! Make sure the ' +
                      'parameter file path is absolute.\nDefaulting to the ' +
                      'parameters in TractLSM/Utils/.')

        return None

    # remove heading '.' in name, which confuses the machine
    if fname[0] == '.':
        fname = fname[1:]

    # dynamic import of module containing class
    fname = importlib.import_module(fname)

    for name, obj in inspect.getmembers(fname):

        if inspect.isclass(obj):
            if cname is None:
                cl = obj()

            elif name == cname:
                cl = obj()

    return cl


def read_csv(fname, drop_units=True):

    """
    Reads csv file with two headers, one for the variables, one for
    their units.

    Arguments:
    ----------
    fname: string
        input filename (with path)

    drop_units: boolean
        if True, drops the second row of the columns, i.e. the units

    Returns:
    --------
    df: pandas dataframe
        dataframe containing all the csv data, units header dropped

    columns: array
        original columns names and units present in csv file

    """

    df = (pd.read_csv(fname, header=[0, 1]).dropna(axis=0, how='all')
          .dropna(axis=1, how='all').squeeze())
    columns = df.columns

    if drop_units:  # drop the units (in second header) for readability
        df.columns = df.columns.droplevel(level=1)

        return df, columns

    else:
        return df


def read_netcdf(fname, var_list=None):

    """
    Retrieves netcdf data & stores it into a dataframe

    Arguments:
    ----------
    fname: string
        input filename (with path)

    var_list: array
        variables to slice from the netcdf

    Returns:
    --------
    df: pandas dataframe
        df containing the data

    """

    ds = xr.open_dataset(fname, autoclose=True)  # access the data
    dates = pd.to_datetime(ds.time.values)  # retrieve dates

    # drop grid
    if var_list is None:
        df = ds.squeeze(dim=('x', 'y'), drop=True).to_dataframe()

    else:
        df = ds[var_list].squeeze(dim=('x', 'y'), drop=True).to_dataframe()

    try:  # reindex
        df['dates'] = dates
        df = df.set_index('dates')

    except ValueError:  # there are multiple z dimensions > 1
        pass

    # free memory
    ds.close()
    del ds
    gc.collect()

    return df


def site_info(site_code, site):

    """
    Access info pertaining to a specific FLUXNET site and retrieve lat,
    lon coordinates in two possible ways.

    Arguments:
    ----------
    site_code: string
        FLUXNET site code

    site: string
        FLUXNET site name

    Returns:
    --------
    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    """

    success = True

    try:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        upper_dir = os.path.dirname(os.path.dirname(script_dir))

        while 'src' in upper_dir:
            upper_dir = os.path.dirname(upper_dir)

        fluxsites_folder = os.path.join(os.path.join(upper_dir, 'input'),
                                        'fluxsites')
        metfile = os.path.join(os.path.join(fluxsites_folder, 'met'),
                               "%sFluxnet_met.nc" % (site))
        ds = xr.open_dataset(metfile)
        lat = ds.latitude.values[0][0]
        lon = ds.longitude.values[0][0]

    except Exception:
        try:
            address = "http://sites.fluxdata.org/"

            try:  # python 3+
                urllib.request.urlopen("%s%s/" % (address, site_code),
                                       timeout=10)

            except ImportError:  # python 2.7
                urllib2.urlopen("%s%s/" % (address, site_code), timeout=10)

            site = pd.read_html("http://sites.fluxdata.org/%s/" %
                                (site_code))[0]

        except Exception:  # request to the servers has timed out
            address = "http://www.fluxdata.org:8080/sitepages/siteInfo.aspx?"

            try:
                try:  # python 3+
                    urllib.request.urlopen("%s%s" % (address, site_code),
                                           timeout=10)

                except ImportError:  # python 2.7
                    urllib2.urlopen("%s%s" % (address, site_code), timeout=10)

                site = pd.read_html("%s%s" % (address, site_code))[0]

            except Exception:
                success = False

        if success:  # the request to the servers worked
            site.index = [e.replace(':', '').lower() for e in site.iloc[:, 0]]
            site.index.name = 'general site information'
            site.drop(site.columns[0], axis=1, inplace=True)
            site.drop(index='general site information', inplace=True)

            lat = float(site.loc['latitude'])
            lon = float(site.loc['longitude'])

    return lat, lon
