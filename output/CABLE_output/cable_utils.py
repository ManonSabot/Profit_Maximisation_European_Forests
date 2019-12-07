#!/usr/bin/env python

"""
Various CABLE utilities

That's all folks.

"""
__author__ = ["Martin De Kauwe", "Manon Sabot"]
__version__ = "1.0 (02.08.2018)"
__email__ = ["mdekauwe@gmail.com", "m.e.b.sabot@gmail.com"]


import os
import sys
import netCDF4
import shutil
import tempfile
import numpy as np
import pandas as pd
import calendar


def adjust_nml_file(fname, replacements):

    """
    Adjust the params/flags in the CABLE namelise file. Note this writes
    over whatever file it is given!

    Parameters:
    ----------
    replacements : dictionary
        dictionary of replacement values.

    """

    f = open(fname, 'r')
    param_str = f.read()
    f.close()
    new_str = replace_keys(param_str, replacements)
    fd, path = tempfile.mkstemp()
    os.write(fd, str.encode(new_str))
    os.close(fd)
    shutil.copy(path, fname)
    os.remove(path)


def replace_keys(text, replacements_dict):

    """
    Function expects to find CABLE namelist file formatted key = value.

    Parameters:
    ----------
    text : string
        input file data.
    replacements_dict : dictionary
        dictionary of replacement values.

    Returns:
    --------
    new_text : string
        input file with replacement values
    """

    lines = text.splitlines()

    for i, row in enumerate(lines):

        # skip blank lines
        if not row.strip():
            continue

        elif "=" not in row:
            lines[i] = row
            continue

        elif not row.startswith("&"):
            key = row.split("=")[0]
            val = row.split("=")[1]
            lines[i] = " ".join((key.rstrip(), "=",
                                 replacements_dict.get(key.strip(),
                                 val.lstrip())))

    return '\n'.join(lines) + '\n'


def add_missing_options_to_nml_file(fname, site, line_start=None):

    """
    Some of the flags we may wish to change are missin from the default
    file so we can't adjust them via this script...add them

    """

    if line_start is None:
        line_start = sum(1 for line in open(fname)) - 1

    f = open(fname, "r")
    contents = f.readlines()
    f.close()

    arg = "   cable_user%GS_SWITCH = 'medlyn'\n"
    contents.insert(line_start, arg)
    line_start += 1

    arg = "   cable_user%GW_MODEL = .FALSE.\n"
    contents.insert(line_start, arg)
    line_start += 1

    arg = "   cable_user%or_evap = .TRUE.\n"
    contents.insert(line_start, arg)
    line_start += 1

    tmp_fname = "tmp_%s.nml" % (site)
    f = open(tmp_fname, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    shutil.move(tmp_fname, fname)


def get_svn_info(here, there):

    """
    Add SVN info and cable namelist file to the output file

    """

    os.chdir(there)
    os.system("svn info > tmp_svn")
    fname = 'tmp_svn'
    fp = open(fname, "r")
    svn = fp.readlines()
    fp.close()
    os.remove(fname)

    url = [i.split(":", 1)[1].strip() for i in svn if i.startswith('URL')]
    rev = [i.split(":", 1)[1].strip() for i in svn if i.startswith('Revision')]
    os.chdir(here)

    return url, rev


def add_attributes_to_output_file(nml_fname, fname, url, rev):

    # Add SVN info to output file
    nc = netCDF4.Dataset(fname, 'r+')
    nc.setncattr('cable_branch', url)
    nc.setncattr('svn_revision_number', rev)

    # Add namelist to output file
    fp = open(nml_fname, "r")
    namelist = fp.readlines()
    fp.close()

    for i, row in enumerate(namelist):

        # skip blank lines
        if not row.strip():
            continue

        # Lines without key = value statement
        elif "=" not in row:
            continue

        # Skip lines that are just comments as these can have "=" too
        elif row.strip().startswith("!"):
            continue

        elif not row.startswith("&"):
            key = str(row.strip().split("=")[0]).rstrip()
            val = str(row.strip().split("=")[1]).rstrip()
            nc.setncattr(key, val)

    nc.close()


def ncdump(nc_fid):

    """
    ncdump outputs dimensions, variables and their attribute
    information.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables

    """

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    nc_dims = [dim for dim in nc_fid.dimensions]

    # Variable information.
    nc_vars = [var for var in nc_fid.variables]

    return nc_attrs, nc_dims, nc_vars


def vpd_to_qair(vpd, tair, press):

    """
    Calculates the near surface specific humidity as given in
    Monteith & Unsworth, 1990, knowing the tair, pressure, and VPD.

    Arguments:
    ----------
    vpd: array or float
        saturation vapour pressure [kPa] at Tair

    tair: array or float
        near surface air temperature [deg K]

    press: array or float
        surface air pressure [kPa]

    Returns:
    --------
    The near surface specific humidity [kg.kg-1]

    """

    # saturation vapor pressure (Tetens eq.)
    T = tair - 273.15
    es = 0.61078 * np.exp(17.27 * T / (T + 237.3))  # kPa

    # relative humidity
    RH = 1. - vpd / es  # not %, 0.-1.
    RH = np.where(RH < 0., 0.01 / 100., RH)
    RH = np.where(RH > 1., 1., RH)

    # Qair, with specific gas constant for dry air and for water vapour
    qair = RH * (286.9 * es) / (461.5 * (press - es))  # kg.kg-1

    return qair


def change_site_level(met_fname, site, fixed=None, lai_fname=None):

    new_met_fname = "%s_tmp.nc" % (site)
    shutil.copyfile(met_fname, new_met_fname)
    nc = netCDF4.Dataset(new_met_fname, 'r+')
    (nc_attrs, nc_dims, nc_vars) = ncdump(nc)

    if fixed is not None:
        lai = fixed

    else:
        df = pd.read_csv(lai_fname, header=[0, 1])[('LAI', '[m2.m-2]')]

        # frequency of met forcing: to match LAI clim
        step = nc['time'][1] - nc['time'][0]
        Nhours = step / (60. * 60.)
        repeat_each = int(24 / Nhours)  # repeat LAI values N times
        lai_clim = (df.loc[df.index.repeat(repeat_each)]
                      .reset_index(drop=True))
        lai = lai_clim.copy()
        lai = lai.iloc[:-repeat_each]  # remove extra day (leap years)

        # get number of years in met forcing, to match with the LAI clim
        Nyears = int(round(nc['time'][nc['time'].shape[0]-1] / (step *
                           repeat_each * 365.25)))
        first_year = int(nc['time'].units.split("since ")[1].split("-")[0])
        all_years = [first_year + e for e in range(Nyears)]

        for year in all_years:
            lai = lai.append(lai_clim.iloc[:-repeat_each], ignore_index=True)

            if calendar.isleap(year):  # this time have extra day in
                lai = lai.append(lai_clim.iloc[[-1] * repeat_each],
                                 ignore_index=True)

        # rm all extra indexes so that it matches the lenght of forcing
        if len(lai) != int(nc['time'].shape[0]):
            lai = lai.iloc[:int(nc['time'].shape[0])]

    LAI = nc.createVariable('LAI', 'f4', ('time', 'y', 'x'))
    LAI.units = "-"
    LAI.long_name = "Leaf Area Index"
    nc.variables['LAI'][:, 0, 0] = lai.values

    cdir = os.path.dirname(os.path.realpath(sys.argv[0]))

    if os.path.isfile(os.path.join(cdir, 'info_sites.csv')):
        df = (pd.read_csv(os.path.join(cdir, 'info_sites.csv'), header=[0, 1])
                .dropna(axis=0, how='all')
                .dropna(axis=1, how='all')
                .squeeze())
        df.columns = [col[0] for col in df.columns]
        df.set_index('Site', inplace=True)

        nc.createVariable('vcmax', 'f4', ('y', 'x'))
        nc.createVariable('ejmax', 'f4', ('y', 'x'))
        nc.createVariable('g1', 'f4', ('y', 'x'))
        nc.createVariable('za', 'f4', ('y', 'x'))
        nc.createVariable('hc', 'f4', ('y', 'x'))

        nc.variables['vcmax'][:] = df.loc[site, 'Vmax25'] * 1e-6
        nc.variables['ejmax'][:] = df.loc[site, 'Vmax25'] * 1.67 * 1e-6
        nc.variables['g1'][:] = df.loc[site, 'g1']
        nc.variables['za'][:] = df.loc[site, 'Zref']  # measurement height
        nc.variables['hc'][:] = df.loc[site, 'Zcan']  # canopy height

    # check that Qair is present in the input
    if 'Qair' not in list(nc.variables):
        Tair = np.asarray(nc.variables['Tair'][:, 0, 0, 0])
        P = np.asarray(nc.variables['PSurf'][:, 0, 0])
        VPD = np.asarray(nc.variables['VPD'][:, 0, 0])

        Qair = nc.createVariable('Qair', 'f8', ('time', 'z', 'y', 'x',))
        Qair.units = "kg/kg"
        Qair.missing_value = -9999.
        Qair.long_name = "Near surface specific humidity"
        Qair.CF_name = "surface_specific_humidity"

        nc.variables['Qair'][:, 0, 0, 0] = vpd_to_qair(VPD, Tair, P)

    nc.close()  # close the new file

    return new_met_fname

