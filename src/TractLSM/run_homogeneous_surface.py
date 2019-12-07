# -*- coding: utf-8 -*-

"""
Run the coupling models between canopy water and carbon fluxes at each
time step given by the forcing file. Soil hydrology is represented by a
simple tipping bucket. The land surface cover is assumed homogeneous.

This file is part of the TractLSM model.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

"""

__title__ = "Run a tractable LSM for a homogeneous surface"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (29.01.2018)"
__email__ = "m.e.b.sabot@gmail.com"


# ======================================================================

# general modules
import collections  # ordered dictionaries
import numpy as np  # array manipulations, math operators
import pandas  # read/write dataframes, csv files

# own modules
from TractLSM import cst  # general constants
from TractLSM.SPAC import wetness, water_potential  # soil modules
from TractLSM.CH2OCoupler import solve_std, set_trans_std  # USO model
from TractLSM.CH2OCoupler import profit_psi  # Profit maximisation

try:  # support functions
    from run_utils import find_model_cases, time_step, write_csv

except (ImportError, ModuleNotFoundError):
    from TractLSM.run_utils import find_model_cases, time_step, write_csv


# ======================================================================

def over_time(df, step, Nsteps, dic, photo, resolution, window, fopt):

    """
    Optimization wrapper at each time step that updates the soil
    moisture and soil water potential for each of the models before
    running them in turn: (i) the Control/USO model (solve_std), (ii)
    the Profit maximisation (profit_psi). None of these are run for
    timesteps when PPFD = 0.

    Arguments:
    ----------
    df: pandas dataframe
        dataframe containing all input data & params

    step: int
        current time step

    Nsteps: int
        total number of steps. This is necessary to know whether unit
        conversion must be based on half-hourly time steps or longer
        time steps!

    dic: dictionary
        initially empty upon input, this dictionary allows to return the
        outputs in a trackable manner. From a time-step to another, it
        also keeps in store the soil moisture and transpiration relative
        to each model, in order to accurately update the soil water
        bucket.

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    resolution: string
        either 'low' (default), 'med', or 'high' to run the optimising
        solver

    fopt: int
        frequency of optimisation

    window: float
        solving window around the last point of optimisation

    Returns:
    --------
    Outputs a tuple of variables depending on the input dic structure.
    When PPFD is zero, a tuple of zero values is returned for A, E, Ci.
    If the models behave in a non-physical manner, zeros are returned
    for A, E, Ci. Overall, the following variables are returned at each
    time step:

        A_std: float
            assimilation rate [umol m-2 s-1] for Ci_std

        A_psi: float
            assimilation rate [umol m-2 s-1] for Ci_psi

        dic['std']['E']: float
            transpiration rate [mmol m-2 s-1] for Tleaf(A_std)

        dic['psi']['E']: float
            transpiration rate [mmol m-2 s-1] for Tleaf(Popt)

        Ci_std: float
            intercellular CO2 concentration [Pa] corresponding to a leaf
            T for which energy balance is met after an iterative solve

        Ci_psi: float
            intercellular CO2 concentration [Pa] corresponding to a leaf
            water potential (P) for which the transpiration cost is
            minimal and the C assimilation gain is maximal

        Rublim_std: string
            'True' if the C assimilation is rubisco limited, 'False'
            otherwise

        Rublim_psi: string
            'True' if the C assimilation is rubisco limited, 'False'
            otherwise

        dic['std']['gs']: float
            stomatal conductance to water [mol m-2 s-1]

        dic['psi']['gs']: float
            stomatal conductance to water [mol m-2 s-1]

        dic['std']['Es']: float
            soil evaporative rate [mmol m-2 s-1]

        dic['psi']['Es']: float
            soil evaporative rate [mmol m-2 s-1]

        dic['std']['sw']: float
            volumetric soil water content [m3 m-3]

        dic['psi']['sw']: float
            volumetric soil water content [m3 m-3]

        dic['std']['Ps']: float
            soil water potential [MPa]

        dic['psi']['Ps']: float
            soil water potential [MPa]

    """

    # parameters & met data
    p = time_step(df, step)

    # tuple of return values
    tpl_return = ()

    # if ProfitMax case number >= 1
    Pcases = [key for key in dic.keys() if 'psi' in key]

    if len(Pcases) >= 1:  # dic to store each ProfitMax case values
        s_cases = [sub.replace('psi', '') for sub in Pcases]
        cProfitMax = {}

    # How many timesteps in a day? (for year, month, hour-based dataset)
    if step == Nsteps - 1:  # last step of time series
        delta = p.hod - df.iloc[step-1, df.columns.get_loc('hod')]

    else:
        if df.iloc[step+1, df.columns.get_loc('hod')] < p.hod:  # end
            delta = df.iloc[step+1, df.columns.get_loc('hod')]

        else:  # during day
            delta = df.iloc[step+1, df.columns.get_loc('hod')] - p.hod

    Dsteps = int(24. / delta)

    # average soil temperature assumed to ~ average air temperature
    if ((step + 1) % Dsteps == 0) and (step != 0):
        Tsoil = (df.iloc[step-(Dsteps-1):step+1, df.columns.get_loc('Tair')]
                   .sum() / float(Dsteps))

    for key in dic.keys():

        # initialise conditions and parameters
        if step == 0:
            try:
                dic[key]['sw'] = p.sw0  # moisture state from spinup

            except AttributeError:
                if p.Ps != p.Psie:  # starting soil moisture from Ps
                    dic[key]['sw'] = water_potential(p, None)

                else:
                    dic[key]['sw'] = p.theta_sat  # saturated soil

            for layer in ['sw0', 'sw1', 'sw2', 'sw3', 'sw4', 'sw5']:

                dic[key][layer] = dic[key]['sw']  # initialise layers

            Tsoil = p.Tair  # simplification
            dic[key]['Tsoil'] = Tsoil  # same Tsoil for every key in dic

            # soil albedo?
            if dic[key]['sw0'] < 0.5 * (p.fc - p.pwp):  # 'dry soil'
                p.albedo_s = p.albedo_ds

            else:  # 'wet soil'
                p.albedo_s = p.albedo_ws

            __, __, __, __, __, __, __, dic[key]['Es'] = \
                wetness(p, Dsteps, dic[key]['sw0'], dic[key]['sw1'],
                        dic[key]['sw2'], dic[key]['sw3'], dic[key]['sw4'],
                        dic[key]['sw5'], 0., 0., Tsoil)

        if ((step + 1) % Dsteps != 0) and (step > 0):
            Tsoil = dic[key]['Tsoil']  # keep same Tsoil through day

        if step > 0:
            dic[key]['Tsoil'] = Tsoil  # same Tsoil for every key in dic

            # soil albedo?
            if dic[key]['sw0'] < 0.5 * (p.fc - p.pwp):  # 'dry soil'
                p.albedo_s = p.albedo_ds

            else:  # 'wet soil'
                p.albedo_s = p.albedo_ws

            dic[key]['sw'], dic[key]['sw0'], dic[key]['sw1'], \
                dic[key]['sw2'], dic[key]['sw3'], dic[key]['sw4'], \
                dic[key]['sw5'], dic[key]['Es'] = wetness(p, Dsteps,
                                                          dic[key]['sw0'],
                                                          dic[key]['sw1'],
                                                          dic[key]['sw2'],
                                                          dic[key]['sw3'],
                                                          dic[key]['sw4'],
                                                          dic[key]['sw5'],
                                                          dic[key]['Es'],
                                                          dic[key]['E'], Tsoil)

        # soil water potential corresponding to the key's soil moisture
        p.Ps = water_potential(p, dic[key]['sw'])
        dic[key]['Ps'] = p.Ps

    # night time
    if (np.isclose(p.PPFD, 0., rtol=cst.zero, atol=cst.zero) or
       np.isclose(0., p.PPFD, rtol=cst.zero, atol=cst.zero)):
        if 'std' in dic.keys():  # standard/Control model
            dic['std']['E'], dic['std']['gs'], A_std, Ci_std, Rublim_std = \
                (0.,) * 5

        if len(Pcases) >= 1:  # ProfitMax

            for icase in range(len(s_cases)):

                PSI = 'psi%s' % (s_cases[icase])
                dic[PSI]['E'], dic[PSI]['gs'], cProfitMax['A_%s' % (PSI)], \
                    cProfitMax['Ci_%s' % (PSI)], \
                    cProfitMax['Rublim_%s' % (PSI)] = (0.,) * 5

    # day time
    else:
        if 'std' in dic.keys():  # standard/Control model
            if dic['std']['sw'] <= p.pwp:  # veg wilted, no solve
                dic['std']['E'], dic['std']['gs'], A_std, Ci_std, \
                    Rublim_std = (0.,) * 5

            else:
                try:
                    p.Ps = dic['std']['Ps']  # use right Ps

                    # soil albedo?
                    if dic['std']['sw0'] < 0.5 * (p.fc - p.pwp):
                        p.albedo_s = p.albedo_ds

                    else:
                        p.albedo_s = p.albedo_ws

                    dic['std']['E'], dic['std']['gs'], A_std, Ci_std, \
                        Rublim_std = solve_std(p, dic['std']['sw'],
                                               photo=photo)

                    # ensure available water is enough
                    try:
                        next_p = time_step(df, step + 1)

                        # soil albedo?
                        if dic['std']['sw0'] < 0.5 * (p.fc - p.pwp):
                            next_p.albedo_s = p.albedo_ds

                        else:
                            next_p.albedo_s = p.albedo_ws

                        next_sw, __, __, __, __, __, __, __ = \
                            wetness(next_p, Dsteps, dic['std']['sw0'],
                                    dic['std']['sw1'], dic['std']['sw2'],
                                    dic['std']['sw3'], dic['std']['sw4'],
                                    dic['std']['sw5'], dic['std']['Es'],
                                    dic['std']['E'], Tsoil)

                        if next_sw < p.pwp:
                            dic['std']['E'], dic['std']['gs'], A_std, Ci_std, \
                                Rublim_std = set_trans_std(p, Dsteps,
                                                           dic['std']['sw'],
                                                           photo=photo)

                    except KeyError:
                        pass

                except (TypeError, IndexError, ValueError):  # no solve
                    dic['std']['E'], dic['std']['gs'], A_std, Ci_std, \
                        Rublim_std = (0.,) * 5

        if len(Pcases) >= 1:  # ProfitMax
            if window is not None:
                window = p.window

            if fopt is None:
                onopt = True
                day_clim = None

            else:
                onopt = p.ontrack  # use condition
                day_clim = None

                # use previous day's climatology to 'aim' for fstom
                if (all(df.iloc[0:Dsteps, df.columns.get_loc('ontrack')]) and
                   (step >= Dsteps)):
                    day_clim = True

                    if onopt:
                        day_clim = (df[df.iloc[step-Dsteps:step,
                                       df.columns.get_loc('PPFD')] >= cst.zero]
                                    .mean())

                        # soil albedo?
                        if dic[PSI]['sw0'] < 0.5 * (p.fc - p.pwp):
                            day_clim.albedo_s = p.albedo_ds

                        else:
                            day_clim.albedo_s = p.albedo_ws

            for icase in range(len(s_cases)):

                PSI = 'psi%s' % (s_cases[icase])
                this_case = int(s_cases[icase])

                try:
                    p.Ps = dic[PSI]['Ps']  # use right Ps

                    # soil albedo?
                    if dic[PSI]['sw0'] < 0.5 * (p.fc - p.pwp):
                        p.albedo_s = p.albedo_ds

                    else:
                        p.albedo_s = p.albedo_ws

                    # get previous optimum
                    isun = df.columns.get_loc('fstom_clim_sun')
                    isha = df.columns.get_loc('fstom_clim_sha')
                    p.fstom_opt_sun = df.iloc[step-1, isun]
                    p.fstom_opt_sha = df.iloc[step-1, isha]

                    if day_clim is not None:
                        if type(day_clim) == bool:  # current optimum
                            p.fstom_opt_sun = df.iloc[step, isun]
                            p.fstom_opt_sha = df.iloc[step, isha]

                        else:  # climatological optimum for future
                            fstom_opt, __, __, __, __, __ = \
                                profit_psi(day_clim, photo=photo,
                                           res=resolution, onopt=onopt,
                                           case=this_case)

                    fstom_opt, dic[PSI]['E'], dic[PSI]['gs'], A_psi, Ci_psi, \
                        Rublim_psi = profit_psi(p, photo=photo, res=resolution,
                                                onopt=onopt, case=this_case)

                    try:
                        df.iloc[step:step+Dsteps, isun] = fstom_opt[0]
                        df.iloc[step:step+Dsteps, isha] = fstom_opt[1]

                    except IndexError:
                        df.iloc[step:, isun] = fstom_opt[0]
                        df.iloc[step:, isha] = fstom_opt[1]

                except (TypeError, IndexError, ValueError):  # no solve
                    dic[PSI]['E'], dic[PSI]['gs'], A_psi, Ci_psi, \
                        Rublim_psi = (0.,) * 5

                cProfitMax['A_%s' % (PSI)] = A_psi
                cProfitMax['Ci_%s' % (PSI)] = Ci_psi
                cProfitMax['Rublim_%s' % (PSI)] = Rublim_psi

    if 'std' in dic.keys():
        tpl_return += (A_std, dic['std']['E'], Ci_std, Rublim_std,
                       dic['std']['gs'], dic['std']['Es'], dic['std']['sw'],
                       dic['std']['Ps'],)

    if len(Pcases) >= 1:

        for icase in range(len(s_cases)):

            PSI = 'psi%s' % (s_cases[icase])
            tpl_return += (cProfitMax['A_%s' % (PSI)], dic[PSI]['E'],
                           cProfitMax['Ci_%s' % (PSI)],
                           cProfitMax['Rublim_%s' % (PSI)], dic[PSI]['gs'],
                           dic[PSI]['Es'], dic[PSI]['sw'], dic[PSI]['Ps'],)

    return tpl_return


def run(fname, df, Nsteps, photo, models=['Control', 'ProfitMax2'],
        resolution=None, fopt=None, window=None):

    """
    Runs the profit maximisation algorithm within a simplified LSM,
    alongsite a Control model which follows traditional photosynthesis
    and transpiration coupling schemes.

    Arguments:
    ----------
    fname: string
        output filename

    df: pandas dataframe
        dataframe containing all input data & params

    Nsteps: int
        total number of time steps over which the models will be run

    photo: string
        either the Farquhar model for photosynthesis, or the Collatz
        model

    models: list of strings
        names of the models to call. Calling all the ProfitMax cases can
        be done via the 'ProfitMax12' string or by listing them
        individually. 'ProfitMax' runs the default case = 2.

    resolution: string
        either 'low' (default), 'med', or 'high' to run the optimising
        solver

    fopt: int
        frequency of optimisation

    window: float
        solving window around the last point of optimisation

    Returns:
    --------
    df2: pandas dataframe
        dataframe of the outputs:
            A_std, A_psi, E_std, E_psi, Ci_std, Ci_psi, Rublim_std,
            Rublim_psi, gs_std, gs_psi, Es_std, Es_psi, sw_std, sw_psi,
            Ps_std, Ps_psi

    """

    # two empty dics, to structure the run setup and retrieve the output
    dic = {}  # appropriately run the models
    output_dic = collections.OrderedDict()  # unpack the output in order

    # sub-dic structures
    subdic = {'sw': None, 'sw0': None, 'sw1': None, 'sw2': None, 'sw3': None,
              'sw4': None, 'sw5': None, 'gs': None, 'Ps': None, 'E': None,
              'Es': None, 'Tsoil': None}  # run

    # for the output dic, the order of things matters!
    subdic2 = collections.OrderedDict([('A', None), ('E', None), ('Ci', None),
                                       ('Rublim', None), ('gs', None),
                                       ('Es', None), ('sw', None),
                                       ('Ps', None)])  # output dic

    # create dictionaries of Nones with the right structures
    if ('Control' in models) or ('Control'.lower() in models):
        dic['std'] = subdic.copy()
        output_dic['std'] = subdic2.copy()

    # are any specific ProfitMax cases specified?
    Pcases = find_model_cases(models, 'ProfitMax')

    # if cases aren't specified, then set to the default ProfitMax case
    if (len(Pcases) < 1) and (('ProfitMax' in models) or
       ('ProfitMax'.lower() in models)):
        dic['psi2'] = subdic.copy()
        output_dic['psi2'] = subdic2.copy()

    # if several ProfitMax cases
    if len(Pcases) >= 1:

        for case in Pcases:

            dic['psi%d' % (case)] = subdic.copy()
            output_dic['psi%d' % (case)] = subdic2.copy()

    # how to run the optimisation?
    if resolution is None:
        resolution = 'low'

    if (fopt is not None) or (window is not None):  # like a mask
        df['ontrack'] = 1
        df['ontrack'].where(df['PPFD'] > cst.zero, 0, inplace=True)
        cumsum = (df['ontrack'].cumsum()
                               .sub(df['ontrack'].cumsum()
                                                 .mask(df['ontrack'] != 0)
                                                 .ffill(), fill_value=0)
                               .astype(int))  # daily cumsum daylight

        if window is not None:  # what solving window over P range?
            df['window'] = None
            df['window'].loc[cumsum > 1] = int(window)  # solving window

        if fopt is not None:  # where to exec the opt or not
            df['ontrack'] = False  # initialise for no fopt
            df['ontrack'].loc[(cumsum - 1) % fopt == 0] = True  # opt

            # N time steps in a day?
            if (df.iloc[(1, df.columns.get_loc('hod'))] <
               df.iloc[(0, df.columns.get_loc('hod'))]):
                delta = df.iloc[(1, df.columns.get_loc('hod'))]

            else:
                delta = (df.iloc[(1, df.columns.get_loc('hod'))]
                         - df.iloc[(0, df.columns.get_loc('hod'))])

            Dsteps = int(24. / delta)

            # opt not set to happen often enough: use avg day forcings
            if float(sum(df['ontrack'])) / float(len(df['ontrack'])) < 0.05:
                df.iloc[0:Dsteps, df.columns.get_loc('ontrack')] = True

    # tracker to know what to aim for if freq or window activated!
    df['fstom_clim_sun'] = 0.7  # plants aim for 0.7 at very beginning
    df['fstom_clim_sha'] = 0.7  # plants aim for 0.7 at very beginning

    # soil albedo changes depending on soil wetness
    df['albedo_s'] = df['albedo_ws'].iloc[0]

    # calculate the attributes that won't change in time
    df['soil_volume'] = df['Zbottom'].iloc[0] * df['ground_area'].iloc[0]
    df['soil_top_volume'] = df['Ztop'].iloc[0] * df['ground_area'].iloc[0]

    # non time-sensitive: last valid value propagated until next valid
    df.fillna(method='ffill', inplace=True)

    # run the model(s) over the range of timesteps / the timeseries
    tpl_out = list(zip(*[over_time(df, step, Nsteps, dic, photo, resolution,
                                   window, fopt) for step in range(Nsteps)]))

    # unpack the output tuple 8 by 8 (A, E, Ci, Rublim, gs, Es, sw, Ps)
    track = 0  # initialize

    for key in output_dic.keys():

        output_dic[key]['A'] = tpl_out[track]
        output_dic[key]['E'] = tpl_out[track + 1]
        output_dic[key]['Ci'] = tpl_out[track + 2]
        output_dic[key]['Rublim'] = tpl_out[track + 3]
        output_dic[key]['gs'] = tpl_out[track + 4]
        output_dic[key]['Es'] = tpl_out[track + 5]
        output_dic[key]['sw'] = tpl_out[track + 6]
        output_dic[key]['Ps'] = tpl_out[track + 7]
        track += 8

    # save the outputs to a csv file and get the corresponding dataframe
    df2 = write_csv(fname, df, output_dic)

    return df2
