#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code that performs the selection of the best parameterisation.

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

The  logic is based on PLUMBER:
    * Best, M. J., Abramowitz, G., Johnson, H. R., Pitman, A. J.,
      Balsamo, G., Boone, A., ... & Ek, M. (2015). The plumbing of land
      surface models: benchmarking model performance. Journal of
      Hydrometeorology, 16(3), 1425-1442.

"""

__title__ = "best performance for calibrating any given parameter"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (08.10.2019)"
__email__ = "m.e.b.sabot@gmail.com"


#=======================================================================

import warnings  # ignore these warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# import general modules
import argparse  # read in the user input
import os  # check for files, paths
import sys  # check for files, paths
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files
from scipy import stats  # compute R coeffs and percentile ranks

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM import conv  # unit converter
from TractLSM.Utils import get_main_dir  # locate files
from TractLSM.Utils import read_csv  # read in data files


#=======================================================================

def main(site, project, parameter=None):

    """
    Main function: Generates 'perf_scores.csv' in the base_dir, which
                   contains information on the performance of different
                   configurations of the model, in the form of
                   statistical metrics and ranks. The logic is based on
                   PLUMBER (Best et al., 2015).

    Arguments:
    ----------
    fname: string
        input filename (with path)

    parameter: string
        None is the default

    Returns:
    --------
    'best.txt' in the output project directory.

    """

    # find all related data
    idata, odata = find_site_data(site, project)

    # where to put the 'best' info files
    rproject = os.path.dirname(odata[0])

    if not os.path.isfile(os.path.join(rproject, 'best.txt')):
        f = open(os.path.join(rproject, 'best.txt'), "w+")

    else:
        f = open(os.path.join(rproject, 'best.txt'), "a+")

    if parameter is not None:  # years which account in the calib
        idata = [e for e in idata if (('2002' in e) or ('2003' in e))]
        odata = [e for e in odata if (('2002' in e) or ('2003' in e))]

    # compare the fluxes to the data, calibrate parameter
    best = performance_scores(idata, odata, parameter)

    fname = [e for e in idata if '%s_' % (best) in e][0]
    df, __ = read_csv(fname)

    f.write('%s: %s\n' % (site, best.split(site)[1]))

    if parameter is None:  # Zroot calibration
        f.write('%s: %s\n' % ('rooting depth', str(df.loc[0, 'Zbottom'])))

    if parameter == 'g1':  # g1 calibration
        f.write('%s: %s\n' % ('g1', str(round(df.loc[0, 'g1'], 2))))

    if parameter == 'fw':  # soil moisture stress calibration
        f.write('%s: %s, %s: %s, %s: %s\n' % ('sens. fw function',
                str(round(df.loc[0, 'sfw'], 2)), 'nudge fc',
                str(round(df.loc[0, 'nfc'], 2)), 'nudge pwp',
                str(round(df.loc[0, 'npwp'], 2))))

    f.close()

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def find_site_data(site, project):

    """
    Finds specific site files in both the input and the output for a
    project, including all model configurations.

    Arguments:
    ----------
    site: string
        individual site name to find

    project: string
        name of the repository which contains the data

    Returns:
    --------
    idata: array
        absolute paths to the input data

    odata: array
        absolute paths to the output data

    """

    idata = []
    odata = []
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(get_main_dir())

    # input data dir path
    iname = os.path.join(basedir, 'input')
    iname = os.path.join(os.path.join(iname, 'projects'), project)

    # output data dir path
    oname = os.path.join(basedir, 'output')
    oname = os.path.join(os.path.join(oname, 'projects'), project)

    for file in os.listdir(iname):

        if (file.endswith('.csv')) and ('actual' in file):
            if site in file:
                fname = os.path.join(iname, file)
                idata += [fname]

    for file in os.listdir(oname):

        if (file.endswith('.csv')) and ('actual' in file):
            if site in file:
                fname = os.path.join(oname, file)
                odata += [fname]

    return idata, odata


def find_root(arr):

    """
    Finds substring common to a whole array of strings.

    Arguments:
    ----------
    arr: array
        list of strings amongst which we want to find common elem.

    Returns:
    --------
    res: string
        the common substring for the whole array of strings

    """

    # determine size of the array
    n = len(arr)

    # take first word from array as reference
    s = arr[0]
    l = len(s)

    res = ''

    for i in range(l):

        for j in range(i + 1, l + 1):

            # all possible substrings of our reference string
            stem = s[i:j]
            k = 1

            for k in range(1, n):

                # check if the generated stem is common to all strings
                if stem not in arr[k]:
                    break

            # current substring present in all strings?
            if (k + 1 == n) and (len(res) < len(stem)):
                res = stem

    return res


def performance_scores(idata, odata, parameter):

    """
    Computes statistical measures for different configurations of the
    model (logic based on PLUMBER; Best et al., 2015).

    Arguments:
    ----------
    idata: array
        absolute paths to the input data

    odata: array
        absolute paths to the output data

    parameter: string
        parameter that is being calibrated for

    Returns:
    --------
    best: array
        string names of all the best parameter calibrations

    """

    to_skip = []
    append = False

    if parameter is None:  # Zroot calibration

        for i in range(len(odata)):  # get all the obs and sims per year

            sim, __ = read_csv(odata[i])
            sim.fillna(0., inplace=True)

            # E/ET must be > 70% in order to consider the calib.
            ratio = (sim['E(std)'].sum() /
                     (sim['E(std)'] + sim['Es(std)']).sum())

            if ratio <= 0.7:
                calib = os.path.split(odata[i])[1].split('_')[0]
                to_skip += ['%s_' % (calib)]

        to_skip = np.unique(to_skip)

    # get all the obs and sims per year
    for i in range(len(idata)):

        calib = os.path.split(odata[i])[1].split('_')[0]

        if '%s_' % (calib) not in to_skip:
            obs, __ = read_csv(idata[i])
            sim, __ = read_csv(odata[i])
            obs.fillna(0., inplace=True)
            sim.fillna(0., inplace=True)

            # on a half-hourly basis
            Dsteps = int(24. / (obs.loc[1, 'hod'] - obs.loc[0, 'hod']))

            # restrict to subsets of dates
            year = int(idata[i].split('_')[-1].split('.csv')[0])
            start = pd.to_datetime(year * 1000 + obs['doy'].iloc[0],
                                   format='%Y%j')

            if obs['doy'].iloc[-1] > obs['doy'].iloc[0]:
                end = pd.to_datetime(year * 1000 + obs['doy'].iloc[-1],
                                     format='%Y%j')

            else:
                end = pd.to_datetime((year + 1) * 1000 + obs['doy'].iloc[-1],
                                     format='%Y%j')

            dates = pd.date_range(start=start, end=end, freq='D')
            __, iunique = np.unique(dates.month, return_index=True)
            unique_months = [dates.month[j] for j in sorted(iunique)]

            # all, restrict dates to April-July (~before dry down)
            dates_2_keep = np.in1d(dates.month, unique_months[3:7])

            # restrict to day time and acceptable LAI
            mask = np.logical_and(obs['PPFD'] > 50., obs['LAI'] > 0.001)
            obs = obs.where(mask).fillna(value=0.)
            sim = sim.where(mask).fillna(value=0.)

            # unit conversion
            obs['Qle'] *= conv.Wpm2_2_mmphlfhr
            sim['Qle'] = ((sim['E(std)'] + sim['Es(std)']) *
                          conv.mmolH2Opm2ps_2_mmphlfhr)

            # daily sums
            obs = obs['Qle'].to_frame()
            sim = sim['Qle'].to_frame()
            obs = obs.groupby(obs.index // Dsteps * Dsteps).sum()
            sim = sim.groupby(sim.index // Dsteps * Dsteps).sum()
            obs.reset_index(inplace=True)
            sim.reset_index(inplace=True)
            obs.drop(columns=['index'], inplace=True)
            sim.drop(columns=['index'], inplace=True)

            # restricted dates
            obs = obs[dates_2_keep]
            sim = sim[dates_2_keep]

            # deal with missing or negative values
            mask = np.logical_and(obs > 0., sim > 0.)
            obs = obs[mask]
            sim = sim[mask]
            obs.fillna(0., inplace=True)
            sim.fillna(0., inplace=True)

            # add the info on which calibration is which
            calib = (odata[i].split(os.sep)[-1].split('_')[0])
            obs['calib'] = calib
            sim['calib'] = calib

            # append all per config in a single df
            if not append:
                obss = obs
                sims = sim
                append = True

            else:
                obss = pd.concat([obss, obs], ignore_index=True)
                sims = pd.concat([sims, sim], ignore_index=True)

    try:
        obss.reset_index(inplace=True)
        sims.reset_index(inplace=True)

        # set the index in the order of increasing param tags!
        root = find_root(np.unique(sims['calib'].values))
        pname = [''.join(ee for i, ee in enumerate(e.split(root)[1]) if not
                 (ee.isdigit() and i != 1)) for e in
                 np.unique(sims['calib'].values)]  # rm param number
        pname = np.unique([e for e in pname if len(e) > 1])[0]
        iparam = []

        for i in range(len(np.unique(sims['calib'].values))):

            if i == 0:
                iparam += [root]

            elif parameter is None:  # Zroot calibration
                iparam += ['%s%s%d' % (root, pname, i)]

            else:
                iparam += ['%s%s%d' % (root, pname, i - 1)]

        # performance across years for a given calibration
        df = pd.DataFrame(index=iparam, columns=['NMSE', 'MAE', 'SD', 'P5',
                                                 'P95'])

        for calib in np.unique(sims['calib'].values):

            subobs = obss.loc[obss['calib'] == calib, 'Qle']
            subsim = sims.loc[sims['calib'] == calib, 'Qle']

            if parameter is None:  # change the normalisation for roots
                nmse = np.mean((subsim - subobs) ** 2. * np.mean(subobs) /
                               np.mean(subsim))

            else:
                nmse = np.mean((subsim - subobs) ** 2. / (np.mean(subsim) *
                               np.mean(subobs)))

            mae = np.mean(np.abs(subsim - subobs))
            sd = np.abs(1. - subsim.std() / np.maximum(1.e-9, subobs.std()))
            p5 = np.abs(np.percentile(subsim, 5) - np.percentile(subobs, 5))
            p95 = np.abs(np.percentile(subsim, 95) - np.percentile(subobs, 95))

            if parameter is not None:
                if np.sum(subobs) < 10.:  # exclude ET < 10. over GS
                    nmse = 2.e3
                    mae = 2.e3
                    sd = 2.e3
                    p5 = 2.e3
                    p95 = 2.e3

            df.loc[calib, 'NMSE'] = nmse
            df.loc[calib, 'MAE'] = mae
            df.loc[calib, 'SD'] = sd
            df.loc[calib, 'P5'] = p5
            df.loc[calib, 'P95'] = p95

        # rank the calibrations
        qranks = np.zeros((len(np.unique(df.index)), len(df.columns)))

        for i in range(len(df.columns)):

            x = np.abs(df.loc[df.index, df.columns[i]].values.flatten())
            qranks[:, i] = np.array([stats.percentileofscore(x, a, 'weak')
                                    if not pd.isna(a) else -1. for a in x])

        qranks = np.mean(qranks, axis=1)
        best = df.index[np.argmin(qranks)]

    except Exception:
        best = [e for e in to_skip if not 'Roots' in e][0]
        best = best.split('_')[0]

    return best


#=======================================================================

if __name__ == "__main__":

    # define the argparse settings to read run set up file
    description = ""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('site', type=str, help='site name')
    parser.add_argument('project', type=str,
                        help='project name (data directory name)')
    parser.add_argument('-p', '--parameter', type=str,
                        help='parameter to calibrate')
    args = parser.parse_args()

    main(args.site, args.project, parameter=args.parameter)
