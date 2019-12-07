# -*- coding: utf-8 -*-

"""
Support functions needed to run the TractLSM model.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "useful ancillary run functions"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (24.03.2019)"
__email__ = "m.e.b.sabot@gmail.com"

# ======================================================================

# general modules
import collections  # ordered dictionaries
import pandas as pd  # read/write dataframes, csv files


# ======================================================================

def find_model_cases(models, model_name):

    """
    Sorts out the model cases called by the user. This is relevant for
    the profit maximisation algorithm.

    Arguments:
    ----------
    models: list of strings
        list of the models called by the user

    model_name: string
        particular model and the associated case(s)

    Returns:
    --------
    A list of ints corresponding to each of the cases called upon.

    """

    cases = [sub.lower() for sub in models if model_name.lower() in
             sub.lower()]

    if len(cases) >= 1:
        cases = [sub.replace(model_name.lower(), '') for sub in cases]
        cases = list(filter(None, cases))  # keep case nums only
        cases = [list(sub) for sub in cases]  # list the cases

        # deal with the cases declared bulk, e.g. 16 for cases 1 to 6
        cases1 = [range(int(sub[0]), int(sub[1]) + 1) for sub in cases if
                  len(sub) > 1]
        cases1 = [case for sub in cases1 for case in sub]  # flat list

        # deal with the cases declared each in turn, i.e. coma separated
        cases2 = [case for sub in cases for case in sub if len(sub) < 2]
        cases2 = [int(sub) for sub in cases2]

        # unique case as an int in increasing order, unique
        cases = sorted(set(cases1 + cases2))

        return cases

    if len(cases) < 1:
        return []


def time_step(df, time_step):

    """
    Accesses each row / time-step of a panda dataframe

    Arguments:
    ----------
    df: pandas dataframe
        dataframe containing all input data & params

    time_step: int
        current time step

    Returns:
    --------
    A copy of the time step's data, i.e. the row (pandas series)

    """

    return df.iloc[time_step].copy()


def write_csv(fname, df, dic):

    """
    Writes a csv output file.

    Arguments:
    ----------
    fname: string
        output filename

    df: pandas dataframe
        dataframe containing all input data & params

    dic: ordered dictionary
        dictionary of the outputs returned by the models, it is read in
        order so that the csv created have the right units and headers
        consistently matching the data

    Returns:
    --------
    df2 : pandas dataframe
        dataframe of the outputs:
            A_std, A_psi(case), E_std, E_psi(case), Ci_std,
            Ci_psi(case), Rublim_std, Rublim_psi(case), Es_std,
            Es_psi(case), sw_std, sw_psi(case), Ps_std, Ps_psi(case),
            gs_std, gs_psi(case)

    Also saves the corresponding csv file, in the output/ folder under
    fname.

    """

    # timeseries for wich the optimisation has occured
    len_series = len(list(list(dic.values())[0].values())[0])

    # declare the lists for the column names, units, and data
    columns = ['doy', 'hod']
    units = ['[-]', '[h]']
    valvars = [list(df['doy'])[:len_series], list(df['hod'])[:len_series]]

    # append to the lists the content of the dictionary
    for subkey in list(dic.values())[0].keys():  # loop on output vars

        for key in dic.keys():  # loop on std, psi

            columns += ['%s(%s)' % (subkey, key)]

            if subkey == 'A':
                units += ['[umol m-2 s-1]']

            if subkey == 'E':
                units += ['[mmol m-2 s-1]']

            if subkey == 'Ci':
                units += ['[Pa]']

            if subkey == 'Rublim':
                units += ['[-]']

            if subkey == 'gs':
                units += ['[mol m-2 s-1]']

            if subkey == 'Es':
                units += ['[mmol m-2 s-1]']

            if subkey == 'sw':
                units += ['[m3 m-3]']

            if subkey == 'Ps':
                units += ['[MPa]']

            valvars += [list(dic[key][subkey])]

    df = (pd.DataFrame(valvars)).T
    df.columns = columns
    these_headers = list(zip(df.columns, units))
    df.columns = pd.MultiIndex.from_tuples(these_headers)

    # write the csv
    df.to_csv(fname, index=False, na_rep='', encoding='utf-8')

    return df
