#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code that performs the selection of the best performing runs. Right now,
it is somehow wired to select for best performance within the repository
with "sample" in its name, and across the repositories with "average"
and "extreme" in their names, but that can easily be tweaked.
Performance is also defined as the ability to model the evaporatite flux
between April and mid-August across all years, but that is also easy to
tweak.

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

__title__ = "best performance of calib. kmax and derived from climate"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (08.10.2018)"
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
import shutil  # copy files on the system
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files
from scipy import stats  # compute R coeffs and percentile ranks

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM import conv  # unit converter
from TractLSM.Utils import read_csv  # read in data files


#=======================================================================

def main(fname):

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

    Returns:
    --------
    'perf_scores.csv' in the base_dir, "best" files in the "best"
    folders.

    """

    if not os.path.join('input', 'projects') in fname:
        raise NameError('The file must be stored in input/projects/...\
                         No performance score will be computed.')
        exit(0)

    else:  # computing performance scores for all the configs
        pd.options.mode.chained_assignment = None
        performance_scores(fname)

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def site_params(df, p):

    """
    Updates a pandas series object used to store the model's parameters
    with the site-specific parameters.

    Arguments:
    ----------
    df: pandas dataframe
        dataframe containing all the data

    p: pandas series
        site parameters

    Returns:
    --------
    p: pandas series
        updated model's parameters (with site info)

    """

    ds = df.iloc[0]

    for index in p.index:

        if index in ds.index:
            if p[index] != ds[index]:
                p[index] = ds[index]

    return p


def performance_scores(fname):

    """
    Computes statistical measures for different configurations of the
    model (logic based on PLUMBER; Best et al., 2015).

    Arguments:
    ----------
    fname: string
        input filename (with path)

    Returns:
    --------
    Generates 'perf_scores.csv' in the base_dir and copies files into
    the "best" folders. Also creates a best.txt file containing all the
    relevant info in the "best" folders.

    """

    get_perf = False

    # find all outputs associated with that site around the base project
    site = fname.split(os.path.sep)[-1].split('_')[0]  # site name?

    # all other output directories with same base dir pattern
    base_dir = os.path.dirname(fname).replace(os.path.join('input',
                                                           'projects'),
                                              os.path.join('output',
                                                           'projects'))
    projects = [d for d in os.listdir(os.path.dirname(base_dir)) if
                base_dir.split('projects%s' % (os.path.sep))[-1] in d]
    projects = [d for d in projects if not 'best' in d]
    all_dirs = [base_dir.replace(base_dir.split('projects%s' %
                (os.path.sep))[-1], d) for d in projects]

    # is the performance score file already present in base_dir?
    fperf = os.path.join(base_dir, 'perf_scores.csv')

    if os.path.isfile(fperf):
        df = (pd.read_csv(fperf, header=[0]).dropna(axis=0, how='all')
                                            .dropna(axis=1, how='all')
                                            .squeeze())
        df.index = df['fname']

    else:
        df = (pd.DataFrame([pd.np.nan, pd.np.nan, pd.np.nan, pd.np.nan,
                            pd.np.nan, pd.np.nan, pd.np.nan, pd.np.nan,
                            pd.np.nan, pd.np.nan, pd.np.nan, pd.np.nan,
                            pd.np.nan, pd.np.nan]).T)
        df.columns = ['Site', 'project', 'year', 'NMSE(ET)', 'MAE(ET)',
                      'SD(ET)', 'P5(ET)', 'P95(ET)', 'NMSE(An)', 'MAE(An)',
                      'SD(An)', 'P5(An)', 'P95(An)', 'fname']
        df.index = df['fname']
        get_perf = True

    # get all the obs and sims per year
    for i in range(len(all_dirs)):

        d = all_dirs[i]
        files = sorted([f for f in os.listdir(d) if '.csv' in f])
        unique_years = np.unique([int(file.split('_')[-1].split('.csv')[0]) for
                                  file in files if site in file])

        # remove years after 1st drought, removing structural changes
        if len(unique_years) > 2:
            unique_years = unique_years[:2]

        # append all per config in a single df
        obss1 = None
        obss2 = None
        sims1 = None
        sims2 = None
        trackk = None
        tracky = 0

        for file in files:

            if 'sample' in projects[i]:  # all per config in a single df
                if obss1 is not None:
                    if trackk != file.split('_')[0]:
                        obss1 = None
                        obss2 = None
                        sims1 = None
                        sims2 = None
                        tracky = 0

            if ((site in file) and ('actual' in file) and
               any([str(y) in file for y in unique_years])):
                trackk = file.split('_')[0]
                tracky += 1

                if os.path.isfile(fperf):

                    try:
                        perf1 = df.loc[df['project'] == projects[i]]
                        nmse1 = perf1.loc[file, ['NMSE(ET)']]
                        mae1 = perf1.loc[file, ['MAE(ET)']]
                        sd1 = perf1.loc[file, ['SD(ET)']]
                        p51 = perf1.loc[file, ['P5(ET)']]
                        p951 = perf1.loc[file, ['P95(ET)']]

                        perf2 = df.loc[df['project'] == projects[i]]
                        nmse2 = perf2.loc[file, ['NMSE(An)']]
                        mae2 = perf2.loc[file, ['MAE(An)']]
                        sd2 = perf2.loc[file, ['SD(An)']]
                        p52 = perf1.loc[file, ['P5(An)']]
                        p952 = perf1.loc[file, ['P95(An)']]

                    except (IndexError, KeyError):
                        get_perf = True

                if get_perf:
                    year = int(file.split('_')[-1].split('.csv')[0])

                    # read in the obs and convert Qle to ET (mm/0.5h)
                    af = fname.replace(fname.split('_')[-1].split('.csv')[0],
                                       str(year))
                    df1, __ = read_csv(af)
                    df1.fillna(value=0., inplace=True)

                    # restrict to subsets of dates
                    start = pd.to_datetime(year * 1000 + df1['doy'].iloc[0],
                                           format='%Y%j')

                    if df1['doy'].iloc[-1] > df1['doy'].iloc[0]:
                        end = pd.to_datetime(year * 1000 + df1['doy'].iloc[-1],
                                             format='%Y%j')

                    else:
                        end = pd.to_datetime((year + 1) * 1000 +
                                             df1['doy'].iloc[-1],
                                             format='%Y%j')

                    dates = pd.date_range(start=start, end=end, freq='D')
                    __, iunique = np.unique(dates.month, return_index=True)
                    unique_months = [dates.month[j] for j in sorted(iunique)]

                    # each year, start in April, go to Nov
                    dates_2_score = np.in1d(dates.month, unique_months[3:-2])

                    # all, restrict dates to April-July (~pre dry down)
                    dates_2_keep = np.in1d(dates[dates_2_score].month,
                                           unique_months[3:7])

                    # on a half-hourly basis
                    Dsteps = int(24. / (df1.loc[1, 'hod'] - df1.loc[0, 'hod']))

                    # unit conversion
                    obs1 = df1['Qle'] * conv.Wpm2_2_mmphlfhr
                    obs2 = df1['GPP'] * conv.umolCpm2ps_2_gCpm2phlfhr

                    # restrict to day time and acceptable LAI
                    mask = np.logical_and(df1['PPFD'] > 50.,
                                          df1['LAI'] > 0.001)
                    obs1 = obs1.where(mask).fillna(value=0.)
                    obs2 = obs2.where(mask).fillna(value=0.)

                    # daily sums
                    obs1 = (obs1.groupby(obs1.index // Dsteps * Dsteps).sum()
                                .iloc[dates_2_score])
                    obs2 = (obs2.groupby(obs2.index // Dsteps * Dsteps).sum()
                                .iloc[dates_2_score])

                    # append all per config in a single df
                    if obss1 is not None:
                        obss1 = pd.concat([obss1, obs1.iloc[dates_2_keep]],
                                          ignore_index=True)

                    else:
                        obss1 = obs1.iloc[dates_2_keep]

                    if obss2 is not None:
                        obss2 = pd.concat([obss2, obs2.iloc[dates_2_keep]],
                                          ignore_index=True)

                    else:
                        obss2 = obs2.iloc[dates_2_keep]

                    # read in the sim and apply the right conversions
                    df2, __ = read_csv(os.path.join(d, file))
                    df2.fillna(value=0., inplace=True)

                    try:
                        sim1 = ((df2['E(psi2)'] + df2['Es(psi2)']) *
                                conv.mmolH2Opm2ps_2_mmphlfhr)
                        sim2 = df2['A(psi2)'] * conv.umolCpm2ps_2_gCpm2phlfhr

                    except KeyError:
                        sim1 = ((df2['E(psi1)'] + df2['Es(psi1)']) *
                                conv.mmolH2Opm2ps_2_mmphlfhr)
                        sim2 = df2['A(psi1)'] * conv.umolCpm2ps_2_gCpm2phlfhr

                    sim1 = (sim1.groupby(sim1.index // Dsteps * Dsteps).sum()
                                .iloc[dates_2_score])
                    sim2 = (sim2.groupby(sim2.index // Dsteps * Dsteps).sum()
                                .iloc[dates_2_score])

                    # append all per config in a single df
                    if sims1 is not None:
                        sims1 = pd.concat([sims1, sim1.iloc[dates_2_keep]],
                                          ignore_index=True)

                    else:
                        sims1 = sim1.iloc[dates_2_keep]

                    if sims2 is not None:
                        sims2 = pd.concat([sims2, sim2.iloc[dates_2_keep]],
                                          ignore_index=True)

                    else:
                        sims2 = sim2.iloc[dates_2_keep]

                    mask = np.logical_and(sim1 > 0., obs1 > 0.)
                    sim1 = sim1[mask]
                    sim2 = sim2[mask]

                    nmse1 = (np.mean((sim1 - obs1[mask]) ** 2.) /
                             (np.mean(sim1) * np.mean(obs1[mask])))
                    mae1 = np.mean(np.abs(sim1 - obs1[mask]))
                    sd1 = np.abs(1. - sim1.std() / obs1[mask].std())
                    p51 = np.abs(np.percentile(sim1, 5) -
                                 np.percentile(obs1[mask], 5))
                    p951 = np.abs(np.percentile(sim1, 95) -
                                  np.percentile(obs1[mask], 95))

                    nmse2 = (np.mean((sim2 - obs2[mask]) ** 2.) /
                             (np.mean(sim2) * np.mean(obs2[mask])))
                    mae2 = np.mean(np.abs(sim2 - obs2[mask]))
                    sd2 = np.abs(1. - sim2.std() / obs2[mask].std())
                    p52 = np.abs(np.percentile(sim2, 5) -
                                 np.percentile(obs2[mask], 5))
                    p952 = np.abs(np.percentile(sim2, 95) -
                                  np.percentile(obs2[mask], 95))

                    df = df.append(pd.Series({'Site': site,
                                              'project': projects[i],
                                              'year': year,
                                              'NMSE(ET)': nmse1,
                                              'MAE(ET)': mae1,
                                              'SD(ET)': sd1,
                                              'P5(ET)': p51,
                                              'P95(ET)': p951,
                                              'NMSE(An)': nmse2,
                                              'MAE(An)': mae2,
                                              'SD(An)': sd2,
                                              'P5(An)': p52,
                                              'P95(An)': p952,
                                              'fname': file},
                                   name=file))

                    if '%s_' % (site) in file:
                        sim1 = ((df2['E(std)'] + df2['Es(std)']) *
                                conv.mmolH2Opm2ps_2_mmphlfhr)
                        sim2 = df2['A(std)'] * conv.umolCpm2ps_2_gCpm2phlfhr

                        sim1 = (sim1.groupby(sim1.index // Dsteps * Dsteps)
                                    .sum().iloc[dates_2_score])
                        sim2 = (sim2.groupby(sim2.index // Dsteps * Dsteps)
                                    .sum().iloc[dates_2_score])

                        mask = np.logical_and(sim1 > 0., obs1 > 0.)
                        sim1 = sim1[mask]
                        sim2 = sim2[mask]

                        nmse1 = (np.mean((sim1 - obs1[mask]) ** 2.) /
                                 (np.mean(sim1) * np.mean(obs1[mask])))
                        mae1 = np.mean(np.abs(sim1 - obs1[mask]))
                        sd1 = np.abs(1. - sim1.std() / obs1[mask].std())
                        p51 = np.abs(np.percentile(sim1, 5) -
                                     np.percentile(obs1[mask], 5))
                        p951 = np.abs(np.percentile(sim1, 95) -
                                      np.percentile(obs1[mask], 95))

                        nmse2 = (np.mean((sim2 - obs2[mask]) ** 2.) /
                                 (np.mean(sim2) * np.mean(obs2[mask])))
                        mae2 = np.mean(np.abs(sim2 - obs2[mask]))
                        sd2 = np.abs(1. - sim2.std() / obs2[mask].std())
                        p52 = np.abs(np.percentile(sim2, 5) -
                                     np.percentile(obs2[mask], 5))
                        p952 = np.abs(np.percentile(sim2, 95) -
                                      np.percentile(obs2[mask], 95))

                        df = df.append(pd.Series({'Site': site,
                                                  'project': 'control',
                                                  'year': year,
                                                  'NMSE(ET)': nmse1,
                                                  'MAE(ET)': mae1,
                                                  'SD(ET)': sd1,
                                                  'P5(ET)': p51,
                                                  'P95(ET)': p951,
                                                  'NMSE(An)': nmse2,
                                                  'MAE(An)': mae2,
                                                  'SD(An)': sd2,
                                                  'P5(An)': p52,
                                                  'P95(An)': p952,
                                                  'fname': file},
                                                 name=file))

            # overall performance across years, wet years, and dry years
            if 'sample' in projects[i]:
                if 'actual' in file:
                    if tracky >= len(unique_years):
                        mask = np.logical_and(sims1 > 0., obss1 > 0.)
                        obss1 = obss1[mask]
                        obss2 = obss2[mask]
                        sims1 = sims1[mask]
                        sims2 = sims2[mask]

                        nmse1 = (np.mean((sims1 - obss1) ** 2.) /
                                 (np.mean(sims1) * np.mean(obss1)))
                        mae1 = np.mean(np.abs(sims1 - obss1))
                        sd1 = np.abs(1. - sims1.std() / obss1.std())
                        p51 = np.abs(np.percentile(sims1, 5) -
                                     np.percentile(obss1, 5))
                        p951 = np.abs(np.percentile(sims1, 95) -
                                      np.percentile(obss1, 95))

                        nmse2 = (np.mean((sims2 - obss2) ** 2.) /
                                 (np.mean(sims2) * np.mean(obss2)))
                        mae2 = np.mean(np.abs(sims2 - obss2))
                        sd2 = np.abs(1. - sims2.std() / obss2.std())
                        p52 = np.abs(np.percentile(sims2, 5) -
                                     np.percentile(obss2, 5))
                        p952 = np.abs(np.percentile(sims2, 95) -
                                      np.percentile(obss2, 95))

                        df = df.append(pd.Series({'Site': site,
                                                  'project': projects[i],
                                                  'year': 'overall',
                                                  'NMSE(ET)': nmse1,
                                                  'MAE(ET)': mae1,
                                                  'SD(ET)': sd1,
                                                  'P5(ET)': p51,
                                                  'P95(ET)': p951,
                                                  'NMSE(An)': nmse2,
                                                  'MAE(An)': mae2,
                                                  'SD(An)': sd2,
                                                  'P5(An)': p52,
                                                  'P95(An)': p952,
                                                  'fname': file.split('_')[0]},
                                                 name=file.split('_')[0]))

        # overall performance across years, wet years, and dry years
        if 'sample' not in projects[i]:
            if tracky >= len(unique_years) * 3:
                mask = np.logical_and(sims1 > 0., obss1 > 0.)
                obss1 = obss1[mask]
                obss2 = obss2[mask]
                sims1 = sims1[mask]
                sims2 = sims2[mask]

                nmse1 = (np.mean((sims1 - obss1) ** 2.) /
                         (np.mean(sims1) * np.mean(obss1)))
                mae1 = np.mean(np.abs(sims1 - obss1))
                sd1 = np.abs(1. - sims1.std() / obss1.std())
                p51 = np.abs(np.percentile(sims1, 5) - np.percentile(obss1, 5))
                p951 = np.abs(np.percentile(sims1, 95) -
                              np.percentile(obss1, 95))

                nmse2 = (np.mean((sims2 - obss2) ** 2.) /
                         (np.mean(sims2) * np.mean(obss2)))
                mae2 = np.mean(np.abs(sims2 - obss2))
                sd2 = np.abs(1. - sims2.std() / obss2.std())
                p52 = np.abs(np.percentile(sims2, 5) - np.percentile(obss2, 5))
                p952 = np.abs(np.percentile(sims2, 95) -
                              np.percentile(obss2, 95))

                df = df.append(pd.Series({'Site': site,
                                          'project': projects[i],
                                          'year': 'overall',
                                          'NMSE(ET)': nmse1,
                                          'MAE(ET)': mae1,
                                          'SD(ET)': sd1,
                                          'P5(ET)': p51,
                                          'P95(ET)': p951,
                                          'NMSE(An)': nmse2,
                                          'MAE(An)': mae2,
                                          'SD(An)': sd2,
                                          'P5(An)': p52,
                                          'P95(An)': p952,
                                          'fname': file.split('_')[0]},
                                         name=file.split('_')[0]))

    # do not keep duplicated rows (when there's more than one drought)
    df = df[~df.duplicated()]

    if get_perf:
        df.to_csv(fperf, index=False, na_rep='', encoding='utf-8')

    # find "best" perf scores, i.e. using quantile ranks
    best_dir = '%s_best_calib' % (base_dir)
    cbest_dir = '%s_best_climate' % (base_dir)

    if not os.path.isdir(best_dir):  # make new dirs if they don't exist
        os.makedirs(best_dir)

    if not os.path.isdir(cbest_dir):
        os.makedirs(cbest_dir)

    if not os.path.isfile(os.path.join(best_dir, 'best.txt')):
        f1 = open(os.path.join(best_dir, 'best.txt'), "w+")

    else:
        f1 = open(os.path.join(best_dir, 'best.txt'), "a+")

    if not os.path.isfile(os.path.join(cbest_dir, 'best.txt')):
        f2 = open(os.path.join(cbest_dir, 'best.txt'), "w+")

    else:
        f2 = open(os.path.join(cbest_dir, 'best.txt'), "a+")

    # go over the overall perf at the site (for ET only)
    df.index = df['Site']
    pscore = df.loc[site]

    # first, do it for the climate trainings
    ranks = np.zeros((2, 5))  # average & extremes, metrics
    qranks = np.zeros((2, 5))  # average & extremes, metricsc
    average = pscore[pscore['project'].str.contains('adjust_average')]
    extreme = pscore[pscore['project'].str.contains('adjust_extreme')]

    try:
        average = average[average['year'] == 'overall'].filter(regex='(ET)')
        extreme = extreme[extreme['year'] == 'overall'].filter(regex='(ET)')
        average.reset_index(drop=True, inplace=True)
        extreme.reset_index(drop=True, inplace=True)

        idx = np.argsort([(np.abs(average.filter(regex='NMSE'))
                           .values[0][0]),
                          (np.abs(extreme.filter(regex='NMSE'))
                           .values[0][0])])
        rank = 1

        for i in idx:

            ranks[i, 0] = rank
            rank += 1

        x = [(np.abs(average.filter(regex='NMSE')).values[0][0]),
             (np.abs(extreme.filter(regex='NMSE')).values[0][0])]
        qranks[:, 0] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort([np.abs(average.filter(regex='MAE')).values[0][0],
                          np.abs(extreme.filter(regex='MAE')).values[0][0]])
        rank = 1

        for i in idx:

            ranks[i, 1] = rank
            rank += 1

        x = [(np.abs(average.filter(regex='MAE')).values[0][0]),
             (np.abs(extreme.filter(regex='MAE')).values[0][0])]
        qranks[:, 1] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort([average.filter(regex='SD').values[0][0],
                          extreme.filter(regex='SD').values[0][0]])
        rank = 1

        for i in idx:

            ranks[i, 2] = rank
            rank += 1

        x = [(np.abs(average.filter(regex='SD')).values[0][0]),
             (np.abs(extreme.filter(regex='SD')).values[0][0])]
        qranks[:, 2] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort([average.filter(regex='P5').values[0][0],
                          extreme.filter(regex='P5').values[0][0]])
        rank = 1

        for i in idx:

            ranks[i, 3] = rank
            rank += 1

        x = [(np.abs(average.filter(regex='P5')).values[0][0]),
             (np.abs(extreme.filter(regex='P5')).values[0][0])]
        qranks[:, 3] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort([average.filter(regex='P95').values[0][0],
                          extreme.filter(regex='P95').values[0][0]])
        rank = 1

        for i in idx:

            ranks[i, 4] = rank
            rank += 1

        x = [(np.abs(average.filter(regex='P95')).values[0][0]),
             (np.abs(extreme.filter(regex='P95')).values[0][0])]
        qranks[:, 4] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

    except KeyError:
        pass

    ranks = np.ma.masked_where(ranks <= 0., ranks)
    qranks = np.ma.masked_where(qranks < 0., qranks)

    # average ranking across fluxes
    rank = np.ma.mean(ranks, axis=1)
    qrank = np.ma.mean(qranks, axis=1)

    # copy best perf score output files to new dirs
    if qrank[0] < qrank[1]:
        pro = 'var_kmax_adjust_average'

    else:
        pro = 'var_kmax_adjust_extreme'

    print(site, pro)

    rproject = os.path.join(os.path.dirname(base_dir), pro)

    for file in os.listdir(rproject):

        if (site in file) and ('actual' in file):
            shutil.copy(os.path.join(rproject, file),
                        os.path.join(cbest_dir, file))

    if get_perf:
        f2.write('%s: %s\n' % (site, pro))
        f2.close()

    # second, do it for the calibration
    pscore_calib = pscore[pscore['project'].str.contains('_sample')]

    pscore_calib['fname'] = pscore_calib['fname'].str.split('_', expand=True)
    pscore_calib.set_index('fname', inplace=True)

    # dims are: kn, metrics
    ranks = np.zeros((len(np.unique(pscore_calib.index)), 5))
    qranks = np.zeros((len(np.unique(pscore_calib.index)), 5))

    try:
        scores = (pscore_calib[pscore_calib['year'] == 'overall']
                  .filter(regex='(ET)'))
        scores.sort_index(inplace=True)

        NMSE = scores.filter(regex='NMSE')
        NMSE.sort_index(inplace=True)
        MAE = scores.filter(regex='MAE')
        MAE.sort_index(inplace=True)
        SD = scores.filter(regex='SD')
        SD.sort_index(inplace=True)
        P5 = scores.filter(regex='P5')
        P5.sort_index(inplace=True)
        P95 = scores.filter(regex='P95')
        P95.sort_index(inplace=True)

        idx = np.argsort(np.abs(NMSE.values.flatten()))
        rank = 1.

        for i in idx:

            ranks[i, 0] = rank
            rank += 0.25

        x = np.abs(NMSE.values.flatten())
        qranks[:, 0] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort(np.abs(MAE.values.flatten()))
        rank = 1.

        for i in idx:

            ranks[i, 1] = rank
            rank += 0.25

        x = np.abs(MAE.values.flatten())
        qranks[:, 1] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort(SD.values.flatten())
        rank = 1.

        for i in idx:

            ranks[i, 2] = rank
            rank += 0.25

        x = np.abs(SD.values.flatten())
        qranks[:, 2] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort(P5.values.flatten())
        rank = 1.

        for i in idx:

            ranks[i, 3] = rank
            rank += 0.25

        x = np.abs(P5.values.flatten())
        qranks[:, 3] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

        idx = np.argsort(P95.values.flatten())
        rank = 1.

        for i in idx:

            ranks[i, 4] = rank
            rank += 0.25

        x = np.abs(P95.values.flatten())
        qranks[:, 4] = np.array([stats.percentileofscore(x, a, 'weak') / 100.
                                 if not pd.isna(a) else -1. for a in x])

    except KeyError:
        pass

    ranks = np.ma.masked_where(ranks <= 0., ranks)
    qranks = np.ma.masked_where(qranks < 0., qranks)

    # average ranking across years, metrics
    rank = np.ma.mean(ranks, axis=1)
    qrank = np.ma.mean(qranks, axis=1)

    best_k = NMSE.index[np.argmin(qrank)].split(site)[1]

    print(site, best_k)

    rproject = os.path.join(os.path.dirname(base_dir), 'var_kmax_sample')

    for file in os.listdir(rproject):

        if ('%s%s_' % (site, best_k) in file) and ('actual' in file):
            shutil.copy(os.path.join(rproject, file),
                        os.path.join(best_dir, file))

    if get_perf:
        f1.write('%s: %s\n' % (site, best_k))
        f1.close()

    return


#=======================================================================

if __name__ == "__main__":

    # define the argparse settings to read run set up file
    description = ""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('fname', type=str, help='input base data file')
    args = parser.parse_args()

    main(args.fname)
