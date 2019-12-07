# -*- coding: utf-8 -*-

"""
Take site-level eddy-covariance met data and generate a forcing file.
Right know it is hardwired to the years yr0 and y1, but this can easily
be modified.
Two options:
- Use yr0 = selected year - 1 as the spinup year;
- Or take the average of all years as the spinup if yr0 is not
  available.

In the cases where flux data is also available, this is added to the
forcing file (but will not force the model), for validation purposes
after the model has finished running. If LAI info is also present, this
is added to the forcing file and is used to run the model.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Asner, G. P., Scurlock, M. O., & Hicke, J. A. (2003). Global Ecology &
  Biogeography, 12, 191–205.
* Monteith, J. L., & Unsworth, M. H. (1990). Principles of environmental
  physics. Arnold. SE, London, UK.

"""

__title__ = "Site-level forcings (met, flux, and LAI info in csv)"
__author__ = ["Manon E. B. Sabot", "Martin G. De Kauwe"]
__version__ = "2.0 (05.03.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# import general modules
import os  # check for files, paths
import sys  # check for version on the system
import csv  # write csv files
from itertools import groupby  # if step-change in climatology
import warnings  # catch warnings
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files
from scipy.interpolate import UnivariateSpline  # interpolate and smooth
from scipy.signal import savgol_filter as savgol  # filter outliers out

# own modules
from TractLSM import conv, cst  # unit converter & general constants

try:
    from general_utils import read_csv, read_netcdf  # read in files
    from calculate_solar_geometry import cos_zenith  # solar geometry
    from modis_lai_lat_lon import main as get_modis_lai  # LAI

except (ImportError, ModuleNotFoundError):
    from TractLSM.Utils.general_utils import read_csv, read_netcdf
    from TractLSM.Utils.calculate_solar_geometry import cos_zenith
    from TractLSM.Utils.modis_lai_lat_lon import main as get_modis_lai


# ======================================================================

def main(fname, year=2003, use_real_years=True):

    """
    Main function: reads fluxnet data in netcdf format (open_file),
                   retrieves the forcing of a specific year
                   (generate_actual_forcing) and either retrieves the
                   data of the year before 'actual' or creates a new
                   forcing based on the average data at a site
                   (generate_spinup_forcing). If the fluxes are present,
                   Qle, NEE, and GPP are added as information to the
                   forcing (however, they will not be used to run the
                   model). When possible, a LAI climatology is added.

    Arguments:
    ----------
    fname: string
        csv forcing filename (with path), preferably stored in the
        input/ folder. The corresponding met (and flux and LAI) data
        must be stored in input/fluxsites/.

    year: int
        selected year

    use_real_years: boolean
        if True real spinup data is pulled, if False it is created based
        on the average at the site.

    Returns:
    --------
    A csv file with the name site_met_and_plant_data_actual.csv or
    site_met_and_plant_data_spinup.csv in the input/(projects/project/)
    folder.

    """

    # separate the folder name from the csv file name
    fname_right = ''.join(fname.split(os.path.sep)[-1:])
    fname_left = str(os.path.sep).join(fname.split(os.path.sep)[:-1])
    fname_left2 = fname_left

    # in case this is a project, remove project name
    if ((not fname_left2.endswith('input%s' % str(os.path.sep))) and
       (not fname_left2.endswith('input'))):
        if fname_left2.endswith(str(os.path.sep)):
            fname_left2 = (str(os.path.sep)
                           .join(fname_left2.split(str(os.path.sep))[:-1]))

        fname_left2 = (str(os.path.sep).join((str(os.path.sep)
                       .join(fname_left2.split(str(os.path.sep))[:-1])
                       .split(str(os.path.sep))[:-1])))

    # site name
    site = fname_right.split('_')[:-1]
    site = site[0]

    # fluxnet sites folder path
    fluxsites_folder = os.path.join(fname_left2, 'fluxsites')
    met_folder = os.path.join(fluxsites_folder, 'met')
    flux_folder = os.path.join(fluxsites_folder, 'flux')
    fnames = [os.path.join(met_folder, '%sFluxnet_met.nc' % (site)),
              os.path.join(flux_folder, '%sFluxnet_flux.nc' % (site))]

    # looking for the met data
    if not os.path.isfile(fnames[0]):
        raise NameError("The met data asked for for %s %d is not present in %s\
                         The forcing cannot be generated."
                        % (site, year, met_folder))
        sys.exit(0)

    # is the flux data present?
    if not os.path.isfile(fnames[1]):
        warnings.warn('The flux data for %s %d is not present in %s\
                       The met forcing only will be generated.'
                      % (site, year, flux_folder))
        df2_yr0 = None  # previous year
        df2_yr1 = None  # current year

    # actual or spinup?
    kind = ''.join(fname_right.split('data_')[-1:])
    kind = ''.join(kind.split('_%d.csv' % (year))[:])

    # read in the met forcing data
    df1 = open_file(fnames[0])
    lat = df1.iloc[(0, df1.columns.get_loc('latitude'))]

    if lat <= -30.:
        D0 = '-06-01'

    else:
        D0 = '-01-01'

    df1_yr1 = df1[(df1.index > '%d%s' % (year, D0)) &
                  (df1.index <= '%d%s' % (year + 1, D0))].copy()

    if kind != 'actual':
        df1_yr0 = df1[(df1.index > '%d%s' % (year - 1, D0)) &
                      (df1.index <= '%d%s' % (year, D0))].copy()

        if (len(df1_yr0) == 0) or (len(df1_yr0) < len(df1_yr1)):
            warnings.warn('There is no met data for the spinup year ' +
                          '(e.g. year - 1) for %s %d.\
                          The spinup forcing will be a doy/hour avg.'
                          % (site, year))
            use_real_years = False

    if len(df1_yr1) == 0:
        raise NameError('There is no met data for the actual year ' +
                        '(e.g. year - 1) for %s %d.\
                         Now exiting the actual run.' % (site, year))
        sys.exit(0)

    # read in the flux data
    if os.path.isfile(fnames[1]):
        df2 = open_file(fnames[1])

        df2_yr0 = df2[(df2.index > '%d%s' % (year - 1, D0)) &
                      (df2.index <= '%d%s' % (year, D0))].copy()
        df2_yr1 = df2[(df2.index > '%d%s' % (year, D0)) &
                      (df2.index <= '%d%s' % (year + 1, D0))].copy()

        if (len(df2_yr0) == 0) or (len(df2_yr0) < len(df2_yr1)):
            warnings.warn('There is no flux data for the spinup year ' +
                          '(e.g. year - 1) for %s %d.\
                          If there is no met data either, it will be an avg.'
                          % (site, year))

        if len(df2_yr1) == 0:
            warnings.warn('There is no flux data for the actual year ' +
                          '(e.g. year - 1) for %s %d.' % (site, year))

    if kind == 'actual':
        ofname = os.path.join(fname_left,
                              '%s_met_and_plant_data_actual_%d.csv' % (site,
                                                                       year))
        generate_actual_forcing(ofname, df1_yr1, df2_yr1)

    if kind != 'actual':
        ofname = os.path.join(fname_left,
                              '%s_met_and_plant_data_spinup_%d.csv' % (site,
                                                                       year))
        generate_spinup_forcing(ofname, df1_yr0, df1, df2_yr0, df2,
                                use_real_years)

    df1 = read_csv(ofname, drop_units=False)  # GPP NaNs in the data?

    if int(df1.GPP.isnull().sum()) > 0:  # other GPP data?
        if kind != 'actual':
            fGPP = os.path.join(flux_folder,
                                '%s.%d.synth.hourly.coreplusquality.csv' %
                                (site, year - 1))

        if kind == 'actual':
            fGPP = os.path.join(flux_folder,
                                '%s.%d.synth.hourly.coreplusquality.csv' %
                                (site, year))

        if not os.path.isfile(fGPP):  # site name a bit different?

            try:
                if kind != 'actual':
                    fGPP = [f for f in os.listdir(flux_folder) if (site in f)
                            and (str(year - 1) in f) and
                            ('synth.hourly.coreplusquality.csv' in f)][0]

                if kind == 'actual':
                    fGPP = [f for f in os.listdir(flux_folder) if (site in f)
                            and (str(year) in f) and
                            ('synth.hourly.coreplusquality.csv' in f)][0]

                fGPP = os.path.join(flux_folder, fGPP)

            except IndexError:
                pass

        if not os.path.isfile(fGPP):  # if year missing, use next year

            try:
                if kind != 'actual':
                    fGPP = [f for f in os.listdir(flux_folder) if (site in f)
                            and (str(year) in f) and
                            ('synth.hourly.coreplusquality.csv' in f)][0]

                if kind == 'actual':
                    fGPP = [f for f in os.listdir(flux_folder) if (site in f)
                            and (str(year + 1) in f) and
                            ('synth.hourly.coreplusquality.csv' in f)][0]

                fGPP = os.path.join(flux_folder, fGPP)

            except IndexError:
                pass

        if os.path.isfile(fGPP):  # use core/original GPP if NaNs in GPP
            df2 = pd.read_csv(fGPP, usecols=['DoY', 'Time', 'GPP_f'],
                              dtype={'DoY': np.int32, 'Time': np.float64,
                                     'GPP_f': np.float64})
            columns = ['doy', 'hod', 'GPP']
            units = ['[-]', '[h]', '[umol m-2 s-1]']
            these_headers = list(zip(columns, units))
            df2.columns = pd.MultiIndex.from_tuples(these_headers)

            # replace flux.nc GPP with core/original GPP if less NaNs
            if int(df2.GPP.isnull().sum()) < int(df1.GPP.isnull().sum()):
                df1.drop(('GPP', '[umol m-2 s-1]'), axis=1, inplace=True)
                df = df1.merge(df2, how='inner', on=[('doy', '[-]'),
                                                     ('hod', '[h]')])
                df.to_csv(ofname, index=False, na_rep='', encoding='utf-8')

            else:
                warnings.warn('There are more NaNs in the core GPP data ' +
                              'than in the GPP from the flux.nc file. The ' +
                              'latter data will be kept for %s %d.'
                              % (site, year))

        else:
            warnings.warn('No core GPP data to replace flux.nc GPP NaNs ' +
                          'with for %s %d' % (site, year))

    # is the LAI present?
    LAI_folder = os.path.join(fluxsites_folder, 'LAI')

    # is there not already a climatology?
    fLAI = os.path.join(LAI_folder, '%s_lai_climatology.csv' % (site))

    if not os.path.isfile(fLAI):  # site name a bit different?

        try:
            fLAI = [f for f in os.listdir(LAI_folder) if (site in f) and
                    ('_lai_climatology.csv' in f)][0]
            fLAI = os.path.join(LAI_folder, fLAI)

        except IndexError:
            pass

    if os.path.isfile(fLAI):  # add LAI climatology
        df1 = read_csv(ofname, drop_units=False)  # met & plant data
        df2 = read_csv(fLAI, drop_units=False)  # LAI yearly climatology

        # if user has stored a climatology (with step change?)
        headers = df2.columns  # store original headers
        df2.columns = [col[0] for col in df2.columns]  # drop units

        # count repeated LAI values
        ct = [sum(1 for _ in group) for _, group in groupby(df2.LAI.values)]

        # is the LAI given a step-wise climatology?
        if sum([plateau >= 14 for plateau in ct]) > 8:
            df2.LAI = step_to_linear(df2.LAI)  # interpolate

        df2.columns = headers  # reattach original headers

        # for every doy in df1, append the corresponding climatology
        df = df1.merge(df2, how='inner', on=[('doy', '[-]')])
        df.to_csv(ofname, index=False, na_rep='', encoding='utf-8')

    elif not os.path.isfile(fLAI):  # build climatology if not available
        fLAI = os.path.join(LAI_folder, '%s_filtered_scaled_Lai_500m.csv' %
                            (site))

        if not os.path.isfile(fLAI):  # slightly different site name?

            try:
                fLAI = [f for f in os.listdir(LAI_folder) if (site in f) and
                        ('_filtered_scaled_Lai_500m.csv' in f)][0]
                fLAI = os.path.join(LAI_folder, fLAI)

            except IndexError:
                pass

        if not os.path.isfile(fLAI):  # LAI from MODIS servers: long!!!

            try:
                warnings.warn('LAI will be retrieved from the MODIS servers')
                get_modis_lai(LAI_folder, fluxsites_folder, site=site)

            except ValueError:
                pass

        if os.path.isfile(fLAI):  # add LAI climatology & year LAI
            df1 = read_csv(ofname, drop_units=False)  # met & plant data
            df2 = LAI_climatology(fLAI)  # year LAI climatology

            # for every doy in df1, append the corresponding climatology
            df = df1.merge(df2, how='inner', on=[('doy', '[-]')])
            df.to_csv(ofname, index=False, na_rep='', encoding='utf-8')

            # save the climatology in the LAI folder, to fetching again
            df2.to_csv(os.path.join(LAI_folder,
                                    '%s_lai_climatology.csv' % (site)),
                       index=False, na_rep='', encoding='utf-8')

    return fluxsites_folder


# ======================================================================

# ~~~ Other functions are defined here ~~~

def qair_to_vpd(qair, tair, press):

    """
    Calculates the saturation vapour pressure at a specific temperature
    tair as given in Monteith & Unsworth, 1990. Then calculates the
    actual air vapour pressure, and finally, the VPD.

    Arguments:
    ----------
    qair: array or float
        near surface specific humidity [kg kg-1]

    tair: array or float
        near surface air temperature [degK]

    press: array or float
        surface air pressure [kPa]

    Returns:
    --------
    The vapour pressure deficit [kPa] at T.

    """

    # saturation vapor pressure (Tetens eq.)
    T = tair - conv.C_2_K  # degC
    es = 0.61078 * np.exp(17.27 * T / (T + 237.3))  # kPa

    # actual vapour pressure
    RH = qair * cst.Rv * (press - es) / (cst.Rd * es)  # not %, 0.-1.
    ea = es * RH  # kPa

    return es - ea


def open_file(fname):

    """
    Performs unit conversions on the data, and a variable conversion
    from near surface specific humidity to vapour pressure deficit if
    met data, saves the flux data as it is.

    Arguments:
    ----------
    fname: string
        csv output filename (with path), preferably stored in the input/
        folder. The corresponding fluxnet data must be stored in
        input/fluxsites/

    Returns:
    --------
    df: pandas dataframe
        df containing the met forcing data or flux data

    """

    df = read_netcdf(fname)  # from netcdf to df

    if 'met.nc' in fname:  # solar zenith angle for two-leaf model
        hod = (np.array(df.index.hour).astype(np.float) +
               np.array(df.index.minute).astype(np.float) / 60.)
        df['coszen'] = cos_zenith(np.array(df.index.dayofyear)
                                    .astype(np.float), hod,
                                  df.iloc[(0, df.columns.get_loc('latitude'))],
                                  df.iloc[(0,
                                           df.columns.get_loc('longitude'))])

        df['PPFD'] = df['SWdown'] * conv.SW_2_PAR  # umol m-2 s-1
        df['PPFD'].where(df['PPFD'] >= 50., 0., inplace=True)  # low PAR
        df['PPFD'].where(90. - np.degrees(np.arccos(df['coszen'])) > 0., 0.,
                         inplace=True)  # the sun isn't up

        try:
            df['Rainf'] *= conv.SEC_2_DAY  # mm d-1

        except KeyError:
            df = df.rename(columns={'Precip': 'Rainf'})
            df['Rainf'] *= conv.SEC_2_DAY  # mm d-1

        try:
            df['PSurf'] *= conv.FROM_MILI  # Pa to kPa

        except KeyError:  # use barometric formula
            df['PSurf'] = (101.325 * (df['Tair'] / (df['Tair'] +
                           cst.Lb * df['elevation'])) ** (cst.g0 * cst.Mair *
                           conv.FROM_MILI / (cst.R * cst.Lb)))  # kPa

        try:
            df['VPD'] *= 0.1  # from hPa to kPa

        except KeyError:
            df['VPD'] = qair_to_vpd(df['Qair'], df['Tair'], df['PSurf'])

        df['VPD'].where(df['VPD'] > 0.05, 0.05, inplace=True)  # all sat
        df['Tair'] -= conv.C_2_K  # degC

    return df


def generate_actual_forcing(ofname, df1_yr1, df2_yr1):

    """
    Saves the unaltered data (actual year) in a csv file with the
    appropriate structure.

    Arguments:
    ----------
    ofname: string
        name of the csv which must be saved containing path (default is
        input/)

    df1_yr1: pandas dataframe
        met data for the selected year

    df2_yr1: pandas dataframe
        flux data for the selected year

    Returns:
    --------
    Saves a csv file with a name of the form
    site_met_and_plant_data_actual.csv in the input/ folder.

    """

    # write csv
    ofp = open(ofname, 'w')
    wr = csv.writer(ofp, delimiter=',', quoting=csv.QUOTE_NONE,
                    escapechar=None, dialect='excel')

    # variables and corresponding units for the two headers
    ovars = ['year', 'doy', 'hod', 'coszen', 'PPFD', 'Tair', 'precip', 'VPD',
             'Patm', 'u']
    ounits = ['[-]', '[-]', '[h]', '[-]', '[umol m-2 s-1]', '[deg C]',
              '[mm d-1]', '[kPa]', '[kPa]', '[m s-1]']

    if df2_yr1 is not None:
        ovars += ['Qle', 'NEE', 'GPP']
        ounits += ['[W m-2]', '[umol m-2 s-1]', '[umol m-2 s-1]']

        try:  # from NaN value to NaN
            df2_yr1['GPP'].where(np.logical_and(df2_yr1['GPP'] > -9999.,
                                 df2_yr1['GPP'] < 9999.), inplace=True)

        except KeyError:  # if no GPP variable at all!
            pass

    wr.writerow([var for i, var in enumerate(ovars)])
    wr.writerow([unit for i, unit in enumerate(ounits)])

    # write data
    hod = 0.5

    for i in range(len(df1_yr1)):

        doy = df1_yr1.index.dayofyear[i]

        if hod == 24.:
            doy -= 1

        if df2_yr1 is not None:

            try:
                wr.writerow([df1_yr1.index.year[i], doy, hod,
                             df1_yr1.coszen[i], df1_yr1.PPFD[i],
                             df1_yr1.Tair[i], df1_yr1.Rainf[i], df1_yr1.VPD[i],
                             df1_yr1.PSurf[i], df1_yr1.Wind[i], df2_yr1.Qle[i],
                             df2_yr1.NEE[i], df2_yr1.GPP[i]])

            except AttributeError:
                wr.writerow([df1_yr1.index.year[i], doy, hod,
                             df1_yr1.coszen[i], df1_yr1.PPFD[i],
                             df1_yr1.Tair[i], df1_yr1.Rainf[i], df1_yr1.VPD[i],
                             df1_yr1.PSurf[i], df1_yr1.Wind[i], df2_yr1.Qle[i],
                             df2_yr1.NEE[i]])

        else:
            wr.writerow([df1_yr1.index.year[i], doy, hod, df1_yr1.coszen[i],
                         df1_yr1.PPFD[i], df1_yr1.Tair[i], df1_yr1.Rainf[i],
                         df1_yr1.VPD[i], df1_yr1.PSurf[i], df1_yr1.Wind[i]])

        hod += 1. / float(len(np.unique(df1_yr1.index.minute)))

        if hod == 24.5:
            hod = 0.5

    ofp.close()


def generate_spinup_forcing(ofname, df1_yr0, df1, df2_yr0, df2,
                            use_real_years):

    """
    Saves the either the previous year to actual or the altered data in
    a csv file with the appropriate structure.

    Arguments:
    ----------
    ofname: string
        name of the csv which must be saved containing path (default is
        input/)

    df1_yr0: pandas dataframe
        met data for the selected year - 1

    df1: pandas dataframe
        met data for all years

    df2_yr0: pandas dataframe
        flux data for the selected year - 1

    df2: pandas dataframe
        flux data for all years

    use_real_years: boolean
        if True real spinup data is pulled, if False it is created based
        on the average at the site.

    Returns:
    --------
    Saves a csv file with a name of the form
    site_met_and_plant_data_spinup.csv in the input/ folder.

    """

    # write csv
    ofp = open(ofname, 'w')
    wr = csv.writer(ofp, delimiter=',', quoting=csv.QUOTE_NONE,
                    escapechar=None, dialect='excel')

    # variables and corresponding units for the two headers
    ovars = ['year', 'doy', 'hod', 'coszen', 'PPFD', 'Tair', 'precip', 'VPD',
             'Patm', 'u']
    ounits = ['[-]', '[-]', '[h]', '[-]', '[umol m-2 s-1]', '[deg C]',
              '[mm d-1]', '[kPa]', '[kPa]', '[m s-1]']

    if ((use_real_years) and (df2_yr0 is not None)) or (df2 is not None):
        ovars += ['Qle', 'NEE', 'GPP']
        ounits += ['[W m-2]', '[umol m-2 s-1]', '[umol m-2 s-1]']

        if df2_yr0 is not None:

            try:  # from NaN value to NaN
                df2_yr0['GPP'].where(np.logical_and(df2_yr0['GPP'] > -9999.,
                                     df2_yr0['GPP'] < 9999.), inplace=True)

            except KeyError:  # if no GPP variable at all!
                pass

        if df2 is not None:

            try:  # from NaN value to NaN
                df2['GPP'].where(np.logical_and(df2['GPP'] > -9999.,
                                 df2['GPP'] < 9999.), inplace=True)

            except KeyError:  # if no GPP variable at all!
                pass

    wr.writerow([var for i, var in enumerate(ovars)])
    wr.writerow([unit for i, unit in enumerate(ounits)])

    # write data
    if use_real_years:
        hod = 0.5

        for i in range(len(df1_yr0)):

            doy = df1_yr0.index.dayofyear[i]

            if hod == 24.:
                doy -= 1

            if df2_yr0 is not None:

                try:
                    wr.writerow([df1_yr0.index.year[i], doy, hod,
                                 df1_yr0.coszen[i], df1_yr0.PPFD[i],
                                 df1_yr0.Tair[i], df1_yr0.Rainf[i],
                                 df1_yr0.VPD[i], df1_yr0.PSurf[i],
                                 df1_yr0.Wind[i], df2_yr0.Qle[i],
                                 df2_yr0.NEE[i], df2_yr0.GPP[i]])

                except AttributeError:
                    wr.writerow([df1_yr0.index.year[i], doy, hod,
                                 df1_yr0.coszen[i], df1_yr0.PPFD[i],
                                 df1_yr0.Tair[i], df1_yr0.Rainf[i],
                                 df1_yr0.VPD[i], df1_yr0.PSurf[i],
                                 df1_yr0.Wind[i], df2_yr0.Qle[i],
                                 df2_yr0.NEE[i]])

            else:
                wr.writerow([df1_yr0.index.year[i], doy, hod,
                             df1_yr0.coszen[i], df1_yr0.PPFD[i],
                             df1_yr0.Tair[i], df1_yr0.Rainf[i], df1_yr0.VPD[i],
                             df1_yr0.PSurf[i], df1_yr0.Wind[i]])

            hod += 1. / float(len(np.unique(df1_yr0.index.minute)))

            if hod == 24.5:
                hod = 0.5

        ofp.close()

    else:  # average value at the site for given doy & hour & half hour

        df3 = df1.groupby([df1.index.dayofyear, df1.index.hour,
                           df1.index.minute]).mean()

        if df2 is not None:
            df4 = df2.groupby([df2.index.dayofyear, df2.index.hour,
                               df2.index.minute]).mean()

        Dsteps = 24 * len(np.unique(df1.index.minute))
        year = int(ofname.split('.')[0].split('_')[-1]) - 1
        ndays = df3.index.levels[0][-1]

        # make sure year starts at the right time
        lat = df1.iloc[(0, df1.columns.get_loc('latitude'))]

        if lat <= -30.:
            D0 = '-06-01'

        else:
            D0 = '-01-01'

        new_index = pd.to_datetime(np.arange(0, ndays * 60 * 60 * Dsteps,
                                             60 * 60), unit='s',
                                   origin=pd.Timestamp('%d%s' % (year, D0)))

        df3.index = df3.index.droplevel(0)
        df3.index = new_index
        df4.index = df4.index.droplevel(0)
        df4.index = new_index

        hod = 0.5

        for i in range(len(df3)):

            doy = df3.index.dayofyear[i]

            if hod == 24.:
                doy -= 1

            if df2 is not None:

                try:
                    wr.writerow([df3.index.year[i], doy, hod, df3.coszen[i],
                                 df3.PPFD[i], df3.Tair[i], df3.Rainf[i],
                                 df3.VPD[i], df3.PSurf[i], df3.Wind[i],
                                 df4.Qle[i], df4.NEE[i], df4.GPP[i]])

                except AttributeError:
                    wr.writerow([df3.index.year[i], doy, hod, df3.coszen[i],
                                 df3.PPFD[i], df3.Tair[i], df3.Rainf[i],
                                 df3.VPD[i], df3.PSurf[i], df3.Wind[i],
                                 df4.Qle[i], df4.NEE[i]])

            else:
                wr.writerow([df3.index.year[i], doy, hod, df3.coszen[i],
                             df3.PPFD[i], df3.Tair[i], df3.Rainf[i],
                             df3.VPD[i], df3.PSurf[i], df3.Wind[i]])

            hod += 24. / float(Dsteps)

            if hod == 24.5:
                hod = 0.5

        ofp.close()

    return


def inflect_pts(y):

    """
    Corrects the artifially created "peaks" and "cavities" that arise
    because of the interpolation of the step changes.

    Arguments:
    ----------
    y: numpy array
        interpolated variable displaying "peaks" and "cavities"

    Returns:
    --------
    y: numpy array
        linearly changing variable, no fake oscillations

    """

    diff = np.zeros(y.shape)

    for i in range(len(y) - 1):

        diff[i + 1] = y[i + 1] - y[i]

    idown = np.where(diff < 0)[0]  # index for negative inflection pts

    for i in range(len(idown)):

        if i < len(idown) - 1:  # smooth the peaks/cavities
            y[idown[i]-1:idown[i+1]] = np.linspace(y[idown[i]-1],
                                                   y[idown[i+1]-1],
                                                   idown[i+1] - idown[i] + 1,
                                                   endpoint=True)

    return y


def step_to_linear(y):

    """
    Interpolates linearly between the different plateaus of a step-
    change function (histogram). This artifially creates "peaks" and
    "cavities" that need further correction. The function inflect_pts is
    used to perform that correction.

    Arguments:
    ----------
    y: array or list
        step-changing variable

    Returns:
    --------
    y: numpy array
        linearly changing variable

    """

    y = np.array(y)
    ind_change = np.where(y[:-1] != y[1:])[0] + 1  # step-change indices

    for i in range(len(ind_change)):  # interpolation between all steps

        if i == 0:
            y[:ind_change[i]] = np.linspace(y[0], y[ind_change[i]],
                                            ind_change[i], endpoint=True)

        if i > 0:
            y[ind_change[i-1]:ind_change[i]] = np.linspace(y[i-1],
                                                           y[ind_change[i]],
                                                           ind_change[i] -
                                                           ind_change[i-1],
                                                           endpoint=True)

        if i == len(ind_change) - 1:
            if ind_change[i] != len(y):
                y[ind_change[i]-2:len(y)] = np.linspace(y[ind_change[i]-2],
                                                        y[len(y)-1], len(y) -
                                                        ind_change[i] + 2,
                                                        endpoint=True)

    y = inflect_pts(y)  # deal with artificial peaks/cavities

    return y


def LAI_climatology(fname):

    """
    Calculates the LAI climatology based on MODIS LAI data. As MODIS
    contains NaNs (excluded) a shut-down period of 4/5 days each year
    (excluded), and the product used is with 4/8 day frequency, the
    climatology is shorter than the actual number of days in a year and
    needs to be extrapolated to 366 days.

    Arguments:
    ----------
    fname: string
        name of the csv file for the MODIS data at 500m height,
        containing path

    Returns:
    --------
    A pandas dataframe containing the multi-year climatology.

    """

    # read-in csv column names
    col_names = ['long_name', 'product', 'mod_date', 'lat_lon', 'proc_date',
                 'band', 'lai']
    df = pd.read_csv(fname, names=col_names, na_values=['F'])

    # clean up the dataframe
    date = df['mod_date'].str[1:]
    df['year'] = date.str[:4]
    df['doy'] = date.str[4:]
    df.index = pd.to_datetime(df.year + df.doy, format='%Y%j')
    df = df.drop(['long_name', 'product', 'mod_date', 'lat_lon', 'proc_date',
                  'band'], axis=1)
    ndays = len(np.unique(df.doy))
    freq = int(round(366. / float(ndays)))  # freq of LAI product used

    data = np.zeros(ndays)
    count = np.zeros(ndays)

    for i in range(len(df)):

        if int(df.iloc[(i, df.columns.get_loc('doy'))]) == 1:
            idx = 0

        else:
            try:
                idx += 1

            except UnboundLocalError:  # series doesn't start on Jan 1st
                idx = (int(int(df.iloc[(i, df.columns.get_loc('doy'))]) / freq)
                       - 1)

        val = df.iloc[(i, df.columns.get_loc('lai'))]

        if not np.isnan(val):
            data[idx] += val  # aggregate LAI
            count[idx] += 1.  # how many values for that day?

    # average across years
    count = np.where(count < 1, 1., count)  # account for leap years
    lai = data / count
    outdays = 366  # number of days in output as long as leap years

    # arrange output timeseries length
    xdays = np.arange(1, len(lai) * freq + 1, freq)
    x_extend = np.hstack((xdays - len(lai) * freq, xdays,
                          xdays + len(lai) * freq))  # three stacks
    y_extend = np.tile(lai, 3)  # repeat LAI to match N stack
    weights = np.tile(count, 3)  # apply the right weights to the data

    # spline the time series, smoothing factor depends on the data
    xdoy = np.arange(1, outdays + 1)  # make it as long as leap years

    if (np.mean(lai) < 2.) and (np.amax(lai) < 2.5):  # low LAIs
        smooth = 25.

    else:
        smooth = 5.
        weights[:] = 1.

    spl = UnivariateSpline(x_extend, y_extend, w=weights, k=5, s=smooth)

    # multiple series out to encourage end = start
    lai_out = np.concatenate([spl(xdoy), spl(xdoy), spl(xdoy)])

    # apply 31-day polynomial running mean to rule outliers out
    lai_smooth = savgol(lai_out, 31, 0)[outdays:2*outdays]
    lai_smooth[lai_smooth < 0.] = cst.LAI_thresh  # neg values to 0.001

    # threshold based on Asner et al., 2003
    lai_smooth[lai_smooth > 50.] = cst.LAI_thresh  # set to 0.001

    # new climatology, with appropriate headers
    new_df = pd.DataFrame([xdoy, lai_smooth]).T
    new_df.columns = ['doy', 'LAI']
    units = ['[-]', '[m2 m-2]']
    these_headers = list(zip(new_df.columns, units))
    new_df.columns = pd.MultiIndex.from_tuples(these_headers)
    new_df.doy = new_df.doy.astype(int)  # to match the met & plant data

    return new_df
