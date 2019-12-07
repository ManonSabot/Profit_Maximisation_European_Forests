#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Quantify the uncertainty from 

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have received
along with the TractLSM.

"""

__title__ = ""
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (24.10.2018)"
__email__ = "m.e.b.sabot@gmail.com"

#==============================================================================

import warnings # ignore these warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# import general modules
import argparse # read in the user input
import os # check for files, paths
import sys # check for paths
import numpy as np # array manipulations, math operators
import pandas as pd # read/write dataframes, csv files

import xarray as xr  # read netcdf


#==============================================================================

def main(folder, sites):

    """
    Main: 

    Arguments:
    ----------
    folder: string
        name of the repository which contains the uncertainty files

    sites: array
        names of the sites to plot

    Returns:
    --------

    """

    Tdiff_all = [0., 0.]
    Qdiff_all = [0., 0.]

    for site in sites:

        fname = os.path.join(os.path.join(
                             os.path.dirname(os.path.realpath(sys.argv[0])),
                             folder), "%s.csv" % (site))

        # open the csv
        df = (pd.read_csv(fname, header=[0]).dropna(axis=0, how='all')
              .dropna(axis=1, how='all').squeeze())
        columns = df.columns

        # open the corresponding netcdf
        fname = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                             "%s_out.nc" % (site))
        ds = xr.open_dataset(fname, autoclose=True)  # access the data
        dates = pd.to_datetime(ds.time.values)  # retrieve dates
        sw = ds['SWdown'].squeeze(dim=('x', 'y'), drop=True).to_dataframe()

        SW_2_PAR = 4.57 * 0.5  # SW (W m-2) to PAR (umol m-2 s-1)

        df = df[-len(dates):]  # remove N spinup cycles
        df['dates'] = dates  # date info
        df['PPFD'] = sw.values * SW_2_PAR  # radiation info
        df.set_index('dates', inplace=True)

        # filter the df by the years we're interested in
        df = df[np.logical_or(np.logical_or(df.index.year == 2002,
                df.index.year == 2003), np.logical_or(df.index.year == 2005,
                df.index.year == 2006))]

        # filter by the months we're interested in (i.e. growing season)
        df = df[np.logical_and(df.index.month > 3, df.index.month < 11)]

        # restrict to day time hours
        df = df[df['PPFD'] > 50.]

        # mean absolute differences, SDs, variances, and max differences
        Tdiff = np.mean(np.abs(df['Tair'] - df['Tcan']))
        Qdiff = np.mean(np.abs(df['Qair'] - df['Qcan']))
        Tsd = np.std(np.abs(df['Tair'] - df['Tcan']))
        Qsd = np.std(np.abs(df['Qair'] - df['Qcan']))
        Tvar = np.var(np.abs(df['Tair'] - df['Tcan']))
        Qvar = np.var(np.abs(df['Qair'] - df['Qcan']))
        Tdiffmax = np.amax(np.abs(df['Tair'] - df['Tcan']))
        Qdiffmax = np.amax(np.abs(df['Qair'] - df['Qcan']))

        Tdiff_all[0] += Tdiff
        Qdiff_all[0] += Qdiff

        print(site)
        print('ABSOLUTE Tdiff: %f, Tsd:%f, Tvar:%f, Tdiffmax: %f'
              % (Tdiff, Tsd, Tvar, Tdiffmax))
        print('ABSOLUTE Qdiff: %f, Qsd:%f, Qvar:%f, Qdiffmax: %f'
              % (Qdiff, Qsd, Qvar, Qdiffmax))

        # mean differences, SDs, and max differences
        Tdiff = np.mean(df['Tair'] - df['Tcan'])
        Qdiff = np.mean(df['Qair'] - df['Qcan'])
        Tsd = np.std(df['Tair'] - df['Tcan'])
        Qsd = np.std(df['Qair'] - df['Qcan'])
        Tvar = np.var(df['Tair'] - df['Tcan'])
        Qvar = np.var(df['Qair'] - df['Qcan'])
        Tdiffmax = np.amax(np.abs(df['Tair'] - df['Tcan']))
        Qdiffmax = np.amax(np.abs(df['Qair'] - df['Qcan']))

        Tdiff_all[1] += Tdiff
        Qdiff_all[1] += Qdiff

        print('Tdiff: %f, Tsd:%f, Tvar:%f, Tdiffmax: %f'
              % (Tdiff, Tsd, Tvar, Tdiffmax))
        print('Qdiff: %f, Qsd:%f, Qvar:%f, Qdiffmax: %f'
              % (Qdiff, Qsd, Qvar, Qdiffmax))


    Tdiff_all[0] /= len(sites)
    Qdiff_all[0] /= len(sites)
    Tdiff_all[1] /= len(sites)
    Qdiff_all[1] /= len(sites)

    print('Overall')
    print('ABSOLUTE Tdiff: %f' % (Tdiff_all[0]))
    print('ABSOLUTE Qdiff: %f' % (Qdiff_all[0]))
    print('Tdiff: %f' % (Tdiff_all[1]))
    print('Qdiff: %f' % (Qdiff_all[1]))

    return


#==============================================================================

if __name__ == "__main__":

    # user input
    sites = ['Hyytiala', 'Soroe', 'Loobos', 'Hesse', 'Parco', 'Puechabon', \
             'Rocca1', 'Rocca2', 'ElSaler1', 'Espirra']

    # define the argparse settings to read run set up file
    description = ""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('folder', type=str,
                        help='folder where the uncertainty files are')
    args = parser.parse_args()

    main(args.folder, sites)
        
