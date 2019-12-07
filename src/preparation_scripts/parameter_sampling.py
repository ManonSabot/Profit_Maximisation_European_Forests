#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code that calibrates the Control, e.g. its rooting depth so as to make
it perform reasonably well.

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "samples the parameter inputs for the Control model"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (20.10.2019)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

import warnings  # ignore these warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# import general modules
import argparse  # read in the user input
import os  # check for files, paths
import sys  # make the TractLSM modules loadable
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.Utils import read_csv  # read in data files


# ======================================================================

def main(fname, sampling='Zroot'):

    """
    Main function: takes a forcing file already meeting the model's
                   requirements and updates some of the Control's
                   parameters.

    Arguments:
    ----------
    fname: string
        input filename (with path), must be stored in the input/ folder

    sampling: string
        parameter to sample. By default, it's the effective rooting
        depth.

    Returns:
    --------
    Saves the new data input files in the same repository they
    originated from.

    """

    # read in input file
    df, columns = read_csv(fname)

    if sampling == 'Zroot':
        Z_samples(df, columns, fname)

    if sampling == 'g1':
        g1_samples(df, columns, fname)

    if sampling == 'fw':
        fws_samples(df, columns, fname)

    return


# ======================================================================

# ~~~ Other functions are defined here ~~~

def Z_samples(df, columns, fname):

    """
    Varies the effective rooting depth between half of CABLE's effective
    rooting depth and 1.5 m.

    Arguments:
    ----------
    df: pandas dataframe
        input data which must be modified

    columns: array
        original columns names and units present in df

    fname: string
        input filename (with path), must be stored in the input/ folder

    Returns:
    --------
    Saves the new data input files in the same repository they
    originated from.

    """

    idx = fname.index('_met_and_plant_data_')

    # CABLE's effective rooting depth per site
    froot = np.asarray([0.08379751, 0.1887792, 0.3332133, 0.3167433,
                       0.07643241, 0.0010342])
    zse = np.asarray([0.022, 0.058, 0.154, 0.409, 1.085, 2.872])
    Zeff = np.sum((df.loc[0, 'fc'] - df.loc[0, 'pwp']) * froot * np.sum(zse))
    Zs = np.arange(0.5 * Zeff, 1.55, 0.05)

    for i in range(len(Zs)):

        df.iloc[0, df.columns.get_loc('Zbottom')] = Zs[i]

        if i == 0:
            fname2 = fname

        else:
            fname2 = '%sRoots%d%s' % (fname[:idx], i, fname[idx:])

        df.columns = columns
        df.to_csv(fname2, index=False, na_rep='', encoding='utf-8')

    return


def g1_samples(df, columns, fname):

    """
    Varies g1 between 0.25 * the original site g1 and 4 * the
    original site g1.

    Arguments:
    ----------
    df: pandas dataframe
        input data which must be modified

    columns: array
        original columns names and units present in df

    fname: string
        input filename (with path), must be stored in the input/ folder

    Returns:
    --------
    Saves the new data input files in the same repository they
    originated from.

    """

    idx = fname.index('_met_and_plant_data_')

    g1_ref = df.iloc[0, df.columns.get_loc('g1')]
    g1_low = g1_ref / 4.
    g1_high = g1_ref * 4.
    g1s = np.linspace(g1_low, g1_ref * 0.95, 12)
    g1s = np.concatenate((g1s, np.linspace(g1_ref * 1.05, g1_high, 12)))

    for i in range(len(g1s)):

        df.iloc[0, df.columns.get_loc('g1')] = g1s[i]
        fname2 = '%sg1WUE%d%s' % (fname[:idx], i, fname[idx:])
        df.columns = columns
        df.to_csv(fname2, index=False, na_rep='', encoding='utf-8')

    return


def fws_samples(df, columns, fname):

    """
    Varies the beta functions by creating combinations of variations
    of the sensitivity of beta (between 0.1 and 10), of variations on
    the plant wilting point (between -0.1 and +0.05 of the original
    value), and of variations on the field capacity (between -0.05 and
    +0.1 of the original value).

    Arguments:
    ----------
    df: pandas dataframe
        input data which must be modified

    columns: array
        original columns names and units present in df

    fname: string
        input filename (with path), must be stored in the input/ folder

    Returns:
    --------
    Saves the new data input files in the same repository they
    originated from.

    """

    idx = fname.index('_met_and_plant_data_')

    # sensitivity of the soil moisture stress function
    sfw_ref = df.iloc[0, df.columns.get_loc('sfw')]
    sfw_low = 0.1
    sfw_high = 10.
    sfws = np.linspace(sfw_low, sfw_ref, 10)
    sfws = np.concatenate((sfws, np.linspace(sfw_ref, sfw_high, 10)[1:]))

    # new field capacity & wilting point must not overlap
    fc_ref = df.iloc[0, df.columns.get_loc('fc')]
    pwp_ref = df.iloc[0, df.columns.get_loc('pwp')]
    nfcs = np.asarray([-0.05, 0., 0.05, 0.1])
    npwps = np.asarray([-0.1, -0.05, 0., 0.05])
    ncombis = np.stack(np.meshgrid(nfcs, npwps), -1).reshape(-1, 2)
    ncombis = np.squeeze(ncombis[np.argwhere(ncombis[:, 0] + fc_ref >
                                             ncombis[:, 1] + pwp_ref)])

    track = 0

    for i in range(len(sfws)):

        for j in range(len(ncombis)):

            if not ((sfws[i] == 1.) and (np.sum(ncombis[j, :]) == 0.)):
                df.iloc[0, df.columns.get_loc('sfw')] = sfws[i]
                df.iloc[0, df.columns.get_loc('nfc')] = ncombis[j, 0]
                df.iloc[0, df.columns.get_loc('npwp')] = ncombis[j, 1]

                track += 1

            fname2 = '%sfws%d%s' % (fname[:idx], track - 1, fname[idx:])
            df.columns = columns
            df.to_csv(fname2, index=False, na_rep='', encoding='utf-8')

    return


# ======================================================================

if __name__ == "__main__":

    # define the argparse settings to read run set up file
    description = ""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('fname', type=str, help='input data file name')
    parser.add_argument('-s', '--sampling', type=str, default='roots',
                        help='type of file sampling/calibration to perform')
    args = parser.parse_args()

    main(args.fname, sampling=args.sampling)
