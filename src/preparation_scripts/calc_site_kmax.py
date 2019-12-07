#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code that retrieves maximum hydraulic conductance, assuming long-term
coordination between the photosynthetic traits, the hydraulic traits,
and the climate. Here, the user has the option to chose which type of
climate their interested in, i.e. "average" or "extreme". The input
files are updated accordingly.

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "calculate site-level values of kmax depending on climate"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (08.10.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

import warnings  # ignore these warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# import general modules
import argparse  # read in the user input
import os  # check for files, paths, version on the system
import sys  # make the TractLSM modules loadable
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM import dparams  # model parameters
from TractLSM.Utils import read_csv  # read in data files
from TractLSM.Utils import cru_climate  # local climate
from TractLSM.Utils.build_final_forcings import update_site_params
from TractLSM.SPAC import net_radiation  # required to calculate kmax
from TractLSM.TraitCoordination import optimal_kmax  # calculate kmax


# ======================================================================

def main(fname, photo='Farquhar', strategy='all', VPD=False, Tair=False,
         Txx=False):

    """
    Main function: takes a forcing file already meeting the model's
                   requirements and updates maximum hydraulic
                   conductance depending on local climate.

    Arguments:
    ----------
    fname: string
        input filename (with path), must be stored in the input/ folder

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    strategy: string
        'high' calculates the kmax associated with P12, 'low' the one
        associated with the water potential before the hydraulic cost
        offsets the net profit, and 'optimal' is associated with maximum
        profit. 'all' retrieves all of the above. 'sample' simply varies
        kmax around the original value

    VPD: boolean
        if True, the VPD associated with either Tair or Txx is used to
        calculate kmax

    Tair: boolean
        if True, the "average" climate air temperature is used to
        calculate kmax

    Txx: boolean
        if True, the "extreme" climate air temperature is used to
        calculate kmax

    Returns:
    --------
    Saves the new data input files in the same repository they
    originated from.

    """

    try:
        idx = fname.index('_met_and_plant_data_')

    except ValueError:
        raise NameError('The input file must be stored in input/...\
                         Exiting now')
        exit(0)

    condition1 = (('H_met_and_plant_data_' not in fname) and
                  ('O_met_and_plant_data_' not in fname) and
                  ('L_met_and_plant_data_' not in fname))

    if condition1:
        if strategy == 'high':
            condition2 = os.path.isfile(fname[:idx] + 'H' + fname[idx:])

        if strategy == 'optimal':
            condition2 = os.path.isfile(fname[:idx] + 'O' + fname[idx:])

        if strategy == 'low':
            condition2 = os.path.isfile(fname[:idx] + 'L' + fname[idx:])

        if strategy == 'all':
            condition2 = (os.path.isfile(fname[:idx] + 'H' + fname[idx:]) and
                          os.path.isfile(fname[:idx] + 'O' + fname[idx:]) and
                          os.path.isfile(fname[:idx] + 'L' + fname[idx:]))

    else:
        condition2 = True

    if strategy == 'sample':
        condition1 = True
        condition2 = False

    if (condition1) and (not condition2):
        if any([VPD, Tair, Txx]):
            VPD, Tair = cru_climate(fname[:idx].split(os.sep)[-1], Tair=Tair,
                                    Txx=Txx, VPD=VPD)

            if pd.isna(VPD):
                VPD = False

            if pd.isna(Tair):
                Tair = False

        if not VPD:
            VPD = 1.

        if not Tair:
            Tair = 25.

        p = class2df(dparams)  # make default param class a pd object

        # replace default parameters by those in the site-specific file
        df, columns = read_csv(fname)
        p = update_site_params(df, p)

        if strategy != 'sample':  # add params missing from default
            p = p.append(pd.Series([Tair, VPD, 1., p.albedo_ws],
                                   index=['Tair', 'VPD', 'scale2can',
                                          'albedo_s']))

            try:  # update site params using composite LAI
                p.LAI = cru_climate(fname[:idx].split(os.sep)[-1], Tair=False,
                                    Txx=False, VPD=False, LAI=True)

            except KeyError:
                pass

            # add missing radiation params (no sensitivity to them)
            p = p.append(pd.Series([1000.], index=['PPFD']))
            p = p.append(pd.Series([net_radiation(p)], index=['Rnet']))

            # get optimal kmax
            sol1, sol2, sol3 = optimal_kmax(p, photo, strategy=strategy)

            # replace kmax in site data, save new file with right name
            if strategy == 'all':
                new = [['H', 'O', 'L'], [sol1[0], sol2[0], sol3[0]]]

            if strategy == 'high':
                new = [['H'], [sol1[0]]]

            if strategy == 'optimal':
                new = [['O'], [sol2[0]]]

            if strategy == 'low':
                new = [['L'], [sol3[0]]]

            for i in range(len(new[0])):

                df.iloc[0, df.columns.get_loc('kmax')] = new[1][i]
                fname2 = fname[:idx] + new[0][i] + fname[idx:]
                df.columns = columns
                df.to_csv(fname2, index=False, na_rep='', encoding='utf-8')

        if strategy == 'sample':
            kmax_ref = df.iloc[0, df.columns.get_loc('kmax')]
            kmax_low = kmax_ref / 4.
            kmax_high = kmax_ref * 2.5
            ks = np.linspace(kmax_low, kmax_ref * 0.95, 12)
            ks = np.concatenate((ks, np.linspace(kmax_ref * 1.05, kmax_high,
                                                 12)))

            for i in range(len(ks)):

                df.iloc[0, df.columns.get_loc('kmax')] = ks[i]
                fname2 = '%sk%d%s' % (fname[:idx], i, fname[idx:])
                df.columns = columns
                df.to_csv(fname2, index=False, na_rep='', encoding='utf-8')

    return


# ======================================================================

# ~~~ Other functions are defined here ~~~

def class2df(cl):

    """
    Converts a class to a pandas series.

    """

    attrs = vars(cl)
    cl = {item[0]: item[1] for item in attrs.items()}

    return pd.Series(cl)


if __name__ == "__main__":

    # define the argparse settings to read run set up file
    description = "Run the kmax optimisation for a specific set of \
                   conditions"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('fname', type=str, help='input data file name')
    parser.add_argument('-b', '--behaviour', type=str, default='all',
                        help='behaviour/strategy for kmax optimisation')
    parser.add_argument('-VPD', '--VPD', action='store_true',
                        help='get long term VPD (kPa)')
    parser.add_argument('-Tair', '--Tair', action='store_true',
                        help='get long term Tair (C)')
    parser.add_argument('-Txx', '--Txx', action='store_true',
                        help='get long term 90th percentile Tair (C)')
    args = parser.parse_args()

    main(args.fname, strategy=args.behaviour, VPD=args.VPD, Tair=args.Tair,
         Txx=args.Txx)
