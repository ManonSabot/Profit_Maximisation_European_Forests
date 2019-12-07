#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code that builds forcing data in the right format depending on various
user input!
This is NOT a working script, but it should be rather helpful in order
to build forcing files in the right format.
To make it work, one must change the file names, options, variables, and
function calls directly in the script.

This file is part of the TractLSM project.

Copyright (c) 2019 Manon E. B. Sabot

Please refer to the terms of the MIT License, which you should have
received along with the TractLSM.

References:
-----------
* Asner, G. P., Scurlock, M. O., & Hicke, J. A. (2003). Global Ecology &
  Biogeography, 12, 191–205.
* Christoffersen, B. O., Gloor, M., Fauset, S., Fyllas, N. M.,
  Galbraith, D. R., Baker, T. R., ... & Sevanto, S. (2016). Linking
  hydraulic traits to tropical forest function in a size-structured and
  trait-driven model (TFS v. 1-Hydro).
* Martin‐StPaul, N., Delzon, S., & Cochard, H. (2017). Plant resistance
  to drought depends on timely stomatal closure. Ecology letters,
  20(11), 1437-1447.
* Monteith, J. L., & Unsworth, M. H. (1990). Principles of environmental
  physics. Arnold. SE, London, UK.

"""

__title__ = "Example of how to prepare forcing data"
__author__ = "Manon E. B. Sabot"
__version__ = "1.0 (22.03.2019)"
__email__ = "m.e.b.sabot@gmail.com"


#=======================================================================

# general modules
import os  # check for files, paths
import sys  # check for files, paths
import numpy as np  # array manipulations, math operators
import pandas as pd  # read/write dataframes, csv files

# first make sure that modules can be loaded from TractLSM
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))

# own modules
from TractLSM.Utils import get_main_dir, retrieve_class, read_csv  # read in data
from TractLSM import conv, cst  # unit converter & general constants
from TractLSM import default_params  # default params
from TractLSM.Utils import cos_zenith  # zenithal cosinus

# those are not always needed
from TractLSM.Utils.cru_climate_lat_lon import extreme  # get extreme CRU cl.
from TractLSM.Utils.cru_climate_lat_lon import growing_season  # avg CRU cl.
from TractLSM.Utils.cru_climate_lat_lon import growing_season_mean
from TractLSM.SPAC import net_radiation  # Rnet?
from TractLSM.TraitCoordination import optimal_kmax  # kmax value?

# local modules, support function, not always needed
from preparation_utils import slope, Px, van_Genuchten_to_Clapp_Hornberger
from preparation_utils import tetens_eq  # VPD
from preparation_utils import CO2_from_NOAA  # monthly average CO2
from preparation_utils import build_LAI_climatology  # LAI


#=======================================================================

def main(fparams, fmet, fflux, fclimate, fLAI, ftraits, comps, add_CO2=True,
         add_LAI=True, year=None, climate='avg'):

    """
    Main function: reads data in csv format, retrieves the met forcing
                   of a specific year or of all years, if possible
                   appends flux data (Qle and GPP) as information to the
                   forcing (they will not beused to force the model
                   runs), LAI, soil parameters, and updates site level
                   species information to generate a holistic forcing
                   file. Finally, the met data of the year before
                   'actual' (or the average met data if that's not
                   available) is used to create forcing data for the
                   spinup.

    Arguments:
    ----------
    fparams: string
        local parameter file name (and path), here we suggest it should
        be a .py file containing a default parameter class, but it
        should be equaly easy to load all parameters from a csv file,
        just see below how to!

    fmet: string
        meteorological forcings file (and path) in csv format

    fflux: string
        flux data file name (and path) in csv format, used to evaluate
        model outputs

    fclimate: string
        background climate data file name (and path) in csv format, used
        to retrieve estimates of maximum hydraulic conductance to then
        force the ProfitMax

    fLAI: string
        LAI data file name (and path) in csv format

    ftraits: string
        trait data file name (and path) in csv format

    comps: array
        different component/species names to read in specific LAI and/or
        trait data if running in multicomponent mode

    add_CO2: boolean
        if True, monthly CO2 data is retrieved from the NOAA database
        and made to match the original dataset's length and days 

    add_LAI: boolean
        if True, the LAI will be retrieved and preprocessed from the
        MODIS dataset for the location's coordinate. This switch will
        automatically be False in case fLAI is provided.
        Warning: retrieving MODIS LAI is long!

    year: int
        selected year

    multicomp: boolean
        if True, loops to create each species/component's forcing file

    climate: string
        if None, there is no acclimation (i.e. kmax at 25 degC and 1
        kPa), if"avg" then acclimates to the mean climate, if "xx" to a
        more extreme climate

    Returns:
    --------
    A csv file with a name of the form
    name?#_met_and_plant_data_actual_year.csv or
    name?_met_and_plant_data_spinup_year.csv in the ddata??? folder,
    with # being a spp/component number or nothing and ??? either
    nothing or a specific acclimation.

    """

    # first go through the files that have been declared
    if fparams is None:
        d = default_params()  # general information

    else:
        d = retrieve_class(fparams)

    # wish to rework some parameters? Here is fine... also e.g. lines 303-305
    # e.g. get the b param for the Clapp-Hornberger model from Weibull params
    m = 0.2
    n = 2.75
    d.bch = van_Genuchten_to_Clapp_Hornberger(d, m, n)

    # adding parameters that don't exist in the standard version, for multicomp
    d.zref = 16.9  # height at which the met measurements are made
    d.height = 1.  # height of the species/component, m
    d.ccov = 0.8  # veg cover, unitless
    d.cbare = 0.2  # soil cover, unitless

    # pre-format data frame
    # try to respect this, otherwise it's hard to say whether things will run
    # properly
    columns = ['Patm', 'u', 'CO2', 'O2', 'zref', 'height', 'Vmax25',
               'gamstar25', 'Tref', 'JV', 'Rlref', 'TRlref', 'Kc25', 'Ko25',
               'alpha', 'c1', 'c2', 'eps_l', 'eps_s', 'albedo_l', 'albedo_ws',
               'albedo_ds', 'ccov', 'cbare', 'tau_l', 'chi_l', 'kn', 'Ev',
               'Ej', 'Egamstar', 'Ec', 'Eo', 'deltaSv', 'deltaSj', 'Hdv',
               'Hdj', 'LAI', 'max_leaf_width', 'canopy_intercept', 'P50',
               'P88', 'ratiocrit', 'kmax', 'g1', 'ground_area', 'Ztop',
               'Zbottom', 'Ps', 'Psie', 'hyds', 'theta_sat', 'fc', 'pwp',
               'bch', 'r_soil']

    df = pd.DataFrame([(d.Patm, d.u, d.CO2, d.O2, d.zref, d.height, d.Vmax25,
                        d.gamstar25, d.Tref, d.JV, d.Rlref, d.TRlref, d.Kc25,
                        d.Ko25, d.alpha, d.c1, d.c2, d.eps_l, d.eps_s,
                        d.albedo_l, d.albedo_ws, d.albedo_ds, d.ccov, d.cbare,
                        d.tau_l, d.chi_l, d.kn, d.Ev, d.Ej, d.Egamstar, d.Ec,
                        d.Eo, d.deltaSv, d.deltaSj, d.Hdv, d.Hdj, d.LAI,
                        d.max_leaf_width, d.canopy_intercept, d.P50, d.P88,
                        d.ratiocrit, d.kmax, d.g1, d.ground_area, d.Ztop,
                        d.Zbottom, d.Ps, d.Psie, d.hyds, d.theta_sat, d.fc,
                        d.pwp, d.bch, d.r_soil)], columns=columns)

    # units in final dataset, this is already in the order variables will go in!
    units = ['[-]', '[-]', '[h]', '[-]', '[umol.m-2.s-1]', '[deg C]',
             '[mm.d-1]', '[kPa]', '[kPa]', '[m s-1]', '[W m-2]',
             '[umol s-1 m-2]', '[m2 m-2]', '[Pa]', '[kPa]', '[m]', '[m]',
             '[umol s-1 m-2]', '[Pa]', '[deg C]', '[-]', '[umol s-1 m-2]',
             '[deg C]', '[Pa]', '[Pa]', '[mol(photon) mol(e-)-1]', '[-]',
             '[-]', '[-]', '[-]', '[-]', '[-]', '[-]', '[-]', '[-]', '[-]',
             '[-]', '[-]', '[J mol-1]', '[J mol-1]', '[J mol-1]', '[J mol-1]',
             '[J mol-1]', '[J mol-1 K-1]', '[J mol-1 K-1]', '[J mol-1]',
             '[J mol-1]', '[m2 m-2]', '[m]', '[-]', '[-MPa]', '[-MPa]', '[-]',
             '[mmol m-2 s-1 MPa-1]', '[kPa0.5]', '[m2]', '[m]', '[m]',
             '[MPa]', '[MPa]', '[m s-1]', '[m3 m-3]', '[m3 m-3]', '[m3 m-3]',
             '[-]', '[-]']

    # it's likely the met data needs to be prepared.... but the function might
    # need to be adjusted
    if not add_CO2:
        df1 = prepare_met_data(fmet, d.lat, d.lon, get_CO2=add_CO2)

    else:
        df1 = prepare_met_data(fmet, d.lat, d.lon)

    # read in the fluxes and append them to df1 if possible
    if fflux is not None:
        df2, __ = read_csv(fflux)

        # only keep LE and GPP, but could easily add NEE, etc.
        df2 = df2[['LE', 'GPP']]
        df2.where(df2 > -9999., 0., inplace=True)  # mask invalid data

        # insert the fluxes in the met dataframe (i.e. df1)
        df1.insert(loc=10, column='Qle', value=df2.LE.values[:len(df1)])
        df1.insert(loc=11, column='GPP', value=df2.GPP_f.values[:len(df1)])

    # Now there are two options: 1. we're in multicomponent mode in which case
    # we'll loop over the following, 2. we're not and looping won't happen

    # first, read in the LAI data, whether it's for one component or several
    if fLAI is not None:
        add_LAI = False  # set to False
        df2 = (pd.read_csv(fLAI, header=[0,1]).dropna(axis=0, how='all')
                                              .dropna(axis=1, how='all')
                                              .squeeze())
        df2.columns = df2.columns.droplevel(level=0)

    if add_LAI: # This is VERY long if the MODIS data is not already present
                # (under the filename)!
        df2 = build_LAI_climatology('filename?', lat, lon)

    # read in the plant traits
    if ftraits is not None:
        df3 = (pd.read_csv(ftraits, header=[0]).dropna(axis=0, how='all')
                                               .dropna(axis=1, how='all')
                                               .squeeze())
        df3.set_index('Parameter', inplace=True)  # we might have to set an index

    # in case we're looking at plant acclimation and a climate data file was
    # given, read it in
    if fclimate is not None:
        df4 = (pd.read_csv(fclimate, header=[0]).dropna(axis=0, how='all')
                                                .dropna(axis=1, how='all')
                                                .squeeze())
        tair = df4['tair']  # extract the air temperature
        vpd = tetens_eq(tair, df4['rh'])  # calculate the vpd from rh, kPa

    # now we're moving towards pre-setting hydraulic conductance acclimation!
    if (climate == 'avg') or (climate == 'xx'):
        if fclimate is not None:
            if climate == 'avg':                 
                Tair = np.mean(tair)
                VPD = np.mean(vpd)

            if climate == 'xx':  # 90th percentile Tair & associated VPD 
                N = int(len(df4) / 10)  # 10% largest values
                VPD = vpd[(tair >= tair.nlargest(N).min())]  # VPD largest
                Tair = tair[(tair >= tair.nlargest(N).min())]  # Tair largest
                VPD = VPD.values  # values above and including p90
                Tair = Tair.values
                VPD = VPD[np.argmin(Tair)]  # VPD associated with p90 of Tair
                Tair = np.amin(Tair)  # actual 90th percentile value

        else:  # try to use the CRU datasets if they're stored
            try:
                cr = os.path.join(os.path.join(get_main_dir(), 'input'),
                                  'fluxsites')  # folder where the CRU data is

                if climate == 'avg':  # use the pre-existing retrieval algorithm
                    months = growing_season(cr, lat, lon)  # avg growing season
                    VPD, Tair = growing_season_mean(cr, lat, lon, months)

                if climate == 'xx':  # 90th percentile of Tmax
                    VPD, Tair = extreme(cr, lat, lon)

            except Exception as e:  # no acclimation
                VPD = 1.
                Tair = 25.

    else:  # no acclimation
        VPD = 1.
        Tair = 25.

    if comps is None:  # no multicomponent
        
        if fLAI is not None: # repeat daily LAI to match length of first df
            LAI = np.repeat(df2['LAI'].values, 48)[:len(df1)]  # assuming hlfhr
            LAI[LAI<cst.LAI_thresh] = cst.LAI_thresh  # Asner et al. (2003)

            if add_CO2: # defines the place at which the LAI is added
                df1.insert(loc=12, column='LAI', value=LAI)  # met + fluxes + LAI

            else:
                df1.insert(loc=11, column='LAI', value=LAI)  # met + fluxes + LAI

        if ftraits is not None: # here insert relevant plant traits (once only)
            # e.g. assuming we have local Vcmax data (make sure you use the 
            # same parameter names as in the default parameter file when you 
            # insert things, note that df is the already pre-formatted dataframe 
            df.Vmax25.iloc[0] = df3.loc[('Vcmax25')]

            # here's an example of how to add P88 if you only have P50 or
            # another second value of % hydraulic conductance loss (e.g. P12 is
            # calculated in the following example)
            df.P50.iloc[0] = df3.loc[('P50')]
            P12 = 12. / slope(df.P50.iloc[0])  # using a predefined relationship
            df.P88.iloc[0] = Px(pd.Series([P12, df.P50.iloc[0]],
                                          index=['P12', 'P50']), 88.)

        # retrieve the maximum hydraulic conductance for given plant traits,
        # VPD, Tair, etc. The value of the PPFD doesn't have any incidence so we
        # arbitrarily set it to 1000.
        p = df.iloc[0]  # default forcings with site info 
        p = p.append(pd.Series([Tair, VPD, 1000., 1., d.albedo_ws, d.doy],
                               index=['Tair', 'VPD', 'PPFD', 'scale2can',
                                      'albedo_s', 'doy']))
        p = p.append(pd.Series([net_radiation(p)], index=['Rnet']))

        # here we do it for the Farquhar model by default (that can be changed
        # to Collatz) and we get both the "high" and optimal strategies. The
        # default is the optimal strategy, but the "high" strategy might be more
        # appropriate for certain ecosystems, try it out!
        solh, sol, __ = optimal_kmax(p, 'Farquhar', strategy='all')
        df.kmax.iloc[0] = sol[0]  # use optimal strategy

        # now build the df containing all the prebuilt met + fluxes + plant info
        for ind in df.columns:

            if ind in df1.columns:
                df = df.drop([ind], axis=1)  # drop what's in df1 from df

        if year is not None:
            df1 = df1[df1['year'] == int(year)]
            df2 = df2[df1['year'] == int(year) - 1]  # prepare the spinup data

            if len(df2) == 0:  # no data the year before
                df2 =  # here use a climatology or something

        if year is None:
            df2 = df2[df1['year'] == int(df1['year'].iloc[0])) - 1]

            if len(df2) == 0:  # no data the year before
                df2 =  # here use a climatology or something

        df1 = pd.concat([df1, df], axis=1)
        df2 = pd.concat([df2, df], axis=1)
        these_headers = list(zip(df1.columns, units))  # include units in headers
        df1.columns = pd.MultiIndex.from_tuples(these_headers)
        df2.columns = pd.MultiIndex.from_tuples(these_headers)

        # which project directory must this be stored in?
        ddata = os.path.dirname(fmet)

        if climate == 'avg':  # acclimation directory?
            ddata += '_acclim_mean'

        if climate == 'xx':  # acclimation directory?
            ddata += '_acclim_extremes'

        if not os.path.isdir(ddata):  # if first time writing in acclimation dir
            os.makedirs(ddata)

        if year is None:  # all years, but a single year MUST be in name anyway 
            df1.to_csv(os.path.join(ddata,
                       'name?_met_and_plant_data_actual_%d.csv' %
                       (int(df1['year'].iloc[0]))),
                       index=False, na_rep='', encoding='utf-8')  # save actual
            df2.to_csv(os.path.join(ddata,
                       'name?_met_and_plant_data_spinup_%d.csv' %
                       (int(df1['year'].iloc[0]))),
                       index=False, na_rep='', encoding='utf-8')  # save spinup

        else:  # specific year in final files
            dfa.to_csv(os.path.join(ddata,
                       'name?_met_and_plant_data_actual_%d.csv' %
                       (int(year))),
                       index=False, na_rep='', encoding='utf-8')
            dfc.to_csv(os.path.join(ddata,
                       'name?_met_and_plant_data_spinup_%d.csv' %
                       (int(year))),
                       index=False, na_rep='', encoding='utf-8')

    else: # multicomponent, we're broadly doing the same as before, but with 
          # a few critical tweaks!

        for i in range(len(comps)):  # loop over each spp/component

            dfa = df1.copy()  # copy the already prebuilt met + fluxes

            if fLAI is not None:  # insert the comp's LAI in copy of df1
                LAI = np.repeat(df2[comps[i]].values, 48)[:len(df1)]
                LAI[LAI<cst.LAI_thresh] = cst.LAI_thresh

                if add_CO2:  # defines the place at which the LAI is added
                    dfa.insert(loc=12, column='LAI', value=LAI)

                else:
                    dfa.insert(loc=11, column='LAI', value=LAI)

            if ftraits is not None:  # insert relevant plant traits              
                df.Vmax25.iloc[0] = df3[comps[i]].loc[('Vcmax25')]

                # you can easily reproduce the same calculations as before
                # (e.g. P88, kmax) on a per component basis here

            # build the df containing all the prebuilt met + fluxes + plant info
            dfb = df.copy()

            for ind in dfb.columns:

                if ind in dfa.columns:
                    dfb = dfb.drop([ind], axis=1)

            if year is not None:
                dfa = dfa[dfa['year'] == int(year)]
                dfc = dfa[dfa['year'] == int(year) - 1]

                if len(dfc) == 0:
                    dfc = # here use a climatology or something

            if year is None:
                dfc = dfa[dfa['year'] == int(dfa['year'].iloc[0]) - 1]

                if len(dfc) == 0:
                    dfc =  # here use a climatology or something

            dfa = pd.concat([dfa, dfb], axis=1)
            dfc = pd.concat([dfc, dfb], axis=1)
            these_headers = list(zip(dfa.columns, units))  # units in headers
            dfa.columns = pd.MultiIndex.from_tuples(these_headers)
            dfc.columns = pd.MultiIndex.from_tuples(these_headers)

            # which project directory must this be stored in?
            ddata = os.path.dirname(fmet)

            if (climate == 'avg') and (i == 0):  # acclimation directory?
                ddata = os.path.dirname(fmet)'_acclim_mean'

            if (climate == 'xx') and (i == 0):  # acclimation directory?
                ddata += '_acclim_extremes'

            if not os.path.isdir(ddata):  # first time writing in acclimation dir
                os.makedirs(ddata)

            if year is None:  # all years in final files
                dfa.to_csv(os.path.join(ddata,
                           'name?%d_met_and_plant_data_actual_%d.csv' %
                           (i + 1, int(dfa['year'].iloc[0]))),
                           index=False, na_rep='', encoding='utf-8')
                dfc.to_csv(os.path.join(ddata,
                           'name?%%d_met_and_plant_data_spinup_%d.csv' %
                           (i + 1, int(dfa['year'].iloc[0]))),
                           index=False, na_rep='', encoding='utf-8')

            else:  # specific year in final files
                dfa.to_csv(os.path.join(ddata,
                           'name?%%d_met_and_plant_data_actual_%d.csv' %
                           (i + 1, int(year))),
                           index=False, na_rep='', encoding='utf-8')
                dfc.to_csv(os.path.join(ddata,
                           'name?%%d_met_and_plant_data_spinup_%d.csv' %
                           (i + 1, int(year))),
                           index=False, na_rep='', encoding='utf-8')

    return


#=======================================================================

# ~~~ Other functions are defined here ~~~

def prepare_met_data(fname, lat, lon, get_CO2=True):

    """
    Prepares all the meteorological and radiative forcings, in the right
    units for the TractLSM. The following is just an example for the
    sake of having units and variables in the right spot!
    N.B.: column = 'name' renames a variable in the dataframe accoring
          to our variable names

    Arguments:
    ----------
    fname: string
        name (path) to the original met data

    lat: float
        latitude [decimal degrees]

    lon: float
        longitude [decimal degrees]

    get_CO2: boolean
        if True, monthly CO2 data is retrieved from the NOAA database
        and made to match the original dataset's length and days 

    Returns:
    --------
    dataframe containing the site-level met data and radiative forcings.

    """

    # read in the original data
    df = (pd.read_csv(fname, header=[0]).dropna(axis=0, how='all')
                                        .dropna(axis=1, how='all').squeeze())

    # deal with timeseries, assuming 'year' is declared in the headers (no cap)
    # and at location 0 as well!
    # and that the timesteps are half-hourly... might need reworking
    track = 0

    for yr in np.unique(df['year']):

        if track == 0:
            doy = np.arange(1, 1 + len(df[np.logical_and(df['year'] == yr,
                            df['hour'] == 0)]) / 2)

        if track > 0:
            doy = np.concatenate([doy, np.arange(1, 1 +
                                       len(df[np.logical_and(df['year'] == yr,
                                                             df['hour'] == 0)])
                                       / 2)])

    doy = np.repeat(doy, 48)  # half-hourly timesteps
    df.insert(loc=1, column='doy', value=doy[:len(df)])  # right format doy in df

    # build the hour of day column, assuming hours and minutes are present
    # importantly, the day must start at hod=0.5 and finish at hod=24. This is 
    # a trick for the outputs to match those of a lot of climate models
    hod = df['hour'] + df['mins'] / 60.
    hod += (hod[1] - hod[0])  # must add a timestep at start to finish at hod=24
    df.insert(loc=2, column='hod', value=hod)

    # calculate the cosinus zenithal (for radiation scattering purposes)
    coszen = cos_zenith(doy.astype(np.float), hod, lat, lon)
    df.insert(loc=3, column='coszen', value=coszen)

    # SW to PPFD/PAR, assuming the incoming shortwave component is named SWin
    PPFD = df['Swin'] * conv.SW_2_PAR  # umol.m-2.s-1
    PPFD.where(PPFD >= 50., 0., inplace=True)  # low PAR excluded
    PPFD.where(90. - np.degrees(np.arccos(df['coszen'])) > 0., 0.,
               inplace=True)  # the sun isn't up, exclude low PAR then
    df.insert(loc=4, column='PPFD', value=PPFD)
    
    # rename the air and have it in loc=5
    df.insert(loc=5, column='Tair', value=df['Temp'])  # in deg C

    # mm per hlfhr to mm per day (the unit the model reads rainfall in)
    # so assuming hlfhr we simply multiply by 24 * 2
    df.insert(loc=6, column='precip', value=df['Rainfall'] * 48.)

    # saturation vapor pressure (Tetens eq.)
    vpd = tetens_eq(df['Temp'], df['RH'])  # kPa, provided RH is in %
    df.insert(loc=7, column='VPD', value=vpd)

    df.insert(loc=8, column='Patm', value=df['Pressure'])

    df.insert(loc=9, column='u', value=df['WindSpeed'])

    # retrieves monthly global CO2 from NOAA database?
    if get_CO2:
        CO2 = CO2_from_NOAA(df)
        df.insert(loc=10, column='CO2', value=CO2)

        # drop all other data in the original file that we don't use!
        # Critically, this needs to happen.
        df.drop(df.columns[11:],axis=1,inplace=True)

    else:
        df.drop(df.columns[10:],axis=1,inplace=True)

    return df


#=======================================================================

if __name__ == "__main__":

    
    # directory name (with path) for where the original data is and where to
    # store the model's forcing files. NB: ddata is slightly modified for the
    # acclimation runs (see lines 349-353)
    basedir = get_main_dir()

    while 'src' in basedir:
        basedir = os.path.dirname(basedir)

    ddata = os.path.join(os.path.join(os.path.join(basedir, 'input'),
                                      'projects'), 'project_name?')


    # all relevant input data file names are declared by the user. They can be
    # left blank or to name? if non existing, as that will not prevent the code
    # from running
    fparams = os.path.join(ddata, 'local_params?.py')  # like in projects/test/ 
    fmet = os.path.join(ddata, 'name?.csv')  # met forcings
    fflux = os.path.join(ddata, 'name?.csv')  # flux data
    fclimate = os.path.join(ddata, 'name?.csv')  # long term climate, for acclim
    fLAI = os.path.join(ddata, 'name?.csv')  # lai
    ftraits = os.path.join(ddata, 'name?.csv')  # trait data

    # species or components name list: can be empty if not running multicomp
    comps = ['spp1', 'spp2', 'spp3']

    # check whether the input files exist, if not None will be declared
    if not os.path.isfile(fparams):
        fparams = None

    if not os.path.isfile(fmet):
        print('Met forcings are really needed to create the input \
               files! Please make sure %s is the right absolute path. '
              % (fmet) + 'The script will now exit.')
        exit(1)

    if not os.path.isfile(fflux):
        fflux = None

    if not os.path.isfile(fclimate):
        fclimate = None

    if not os.path.isfile(fLAI):
        fLAI = None

    if not os.path.isfile(ftraits):
        ftraits = None

    if len(comps) == 0:
        comps = None

    main(fparams, fmet, fflux, fclimate, fLAI, ftraits, comps, add_CO2=True,
         year=2011, climate=None)  # no acclimation
    main(fparams, fmet, fflux, fclimate, fLAI, ftraits, comps, add_CO2=True,
         year=2011, climate=None)  # "mean climate" acclim
    main(fparams, fmet, fflux, fclimate, fLAI, ftraits, comps, add_CO2=True,
         year=2011, climate=None)  # "extreme climate" acclim

