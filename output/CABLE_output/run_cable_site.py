#!/usr/bin/env python

"""
Run CABLE either for a single site, a subset, or all the flux sites
pointed to in the met directory

- Only intended for biophysics
- Set mpi = True if doing a number of flux sites

That's all folks.

"""
__author__ = ["Martin De Kauwe", "Manon Sabot"]
__version__ = "1.0 (02.08.2018)"
__email__ = ["mdekauwe@gmail.com", "m.e.b.sabot@gmail.com"]

import os
import glob
import shutil
import multiprocessing as mp
import numpy as np

from cable_utils import adjust_nml_file
from cable_utils import get_svn_info
from cable_utils import change_site_level
from cable_utils import add_missing_options_to_nml_file
from cable_utils import add_attributes_to_output_file


class RunCable(object):

    def __init__(self, met_dir, log_dir, output_dir, restart_dir, aux_dir,
                 namelist_dir, nml_fname, veg_fname, soil_fname, grid_fname,
                 phen_fname, cnpbio_fname, lai_dir, fixed_lai, co2_conc,
                 met_subset, cable_src, cable_exe, mpi, verbose):

        self.met_dir = met_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.aux_dir = aux_dir
        self.namelist_dir = namelist_dir
        self.nml_fname = nml_fname
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")
        self.veg_fname = os.path.join(self.biogeophys_dir, veg_fname)
        self.soil_fname = os.path.join(self.biogeophys_dir, soil_fname)
        self.grid_fname = os.path.join(self.grid_dir, grid_fname)
        self.phen_fname = os.path.join(self.biogeochem_dir, phen_fname)
        self.cnpbio_fname = os.path.join(self.biogeochem_dir, cnpbio_fname)
        self.co2_conc = co2_conc
        self.met_subset = met_subset
        self.cable_src = cable_src
        self.cable_exe = cable_exe
        self.cable_exe = os.path.join(cable_src, "offline/%s" % (cable_exe))
        self.verbose = verbose
        self.mpi = mpi
        self.lai_dir = lai_dir
        self.fixed_lai = fixed_lai

    def main(self, num_cores=None):

        (met_files, url, rev) = self.initialise_stuff()

        # Setup multi-processor jobs
        if self.mpi:
            if num_cores is None:  # use them all!
                num_cores = mp.cpu_count()

            chunk_size = int(np.ceil(len(met_files) / float(num_cores)))
            mp.Pool(processes=num_cores)
            processes = []

            for i in range(num_cores):

                start = chunk_size * i
                end = chunk_size * (i + 1)

                if end > len(met_files):
                    end = len(met_files)

                # setup a list of processes that we want to run
                p = mp.Process(target=self.worker,
                               args=(met_files[start:end], url, rev, ))
                processes.append(p)

            # Run processes
            for p in processes:

                p.start()

            exit_codes = [job.join() for job in processes]
            self.tidy_after_run(met_files, url, rev)

        else:
            self.worker(met_files, url, rev)
            self.tidy_after_run(met_files, url, rev)

    def worker(self, met_files, url, rev):

        for fname in met_files:

            site = os.path.basename(fname).split("Fluxnet")[0]
            print(site)

            base_nml_fn = os.path.join(self.grid_dir, "%s" % (self.nml_fname))
            nml_fname = "cable_%s.nml" % (site)
            shutil.copy(base_nml_fn, nml_fname)
            add_missing_options_to_nml_file(nml_fname, site)

            (out_fname, out_log_fname) = self.clean_up_old_files(site)

            if self.fixed_lai is not None or self.lai_dir is not None:
                if self.lai_dir is not None:
                    lai_fname = os.path.join(self.lai_dir,
                                             "%s_lai_climatology.csv" % (site))

                fname = change_site_level(fname, site, fixed=self.fixed_lai,
                                          lai_fname=lai_fname)

            replace_dict = {
                            "filename%met": "'%s'" % (fname),
                            "filename%out": "'%s'" % (out_fname),
                            "filename%log": "'%s'" % (out_log_fname),
                            "filename%restart_out": "' '",
                            "filename%type": "'%s'" % (self.grid_fname),
                            "filename%veg": "'%s'" % (self.veg_fname),
                            "filename%soil": "'%s'" % (self.soil_fname),
                            "output%restart": ".FALSE.",
                            "fixedCO2": "%.2f" % (self.co2_conc),
                            "casafile%phen": "'%s'" % (self.phen_fname),
                            "casafile%cnpbio": "'%s'" % (self.cnpbio_fname),
                            "cable_user%FWSOIL_SWITCH": "'Haverd2013'",
                            "cable_user%GS_SWITCH": "'medlyn'",
                            "cable_user%or_evap": ".TRUE.",
            }

            adjust_nml_file(nml_fname, replace_dict)

            self.run_me(nml_fname)

    def initialise_stuff(self):

        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        if not os.path.exists(self.namelist_dir):
            os.makedirs(self.namelist_dir)

        # Run all the met files in the directory
        if len(met_subset) == 0:
            met_files = glob.glob(os.path.join(self.met_dir, "*.nc"))

        else:
            met_files = [os.path.join(self.met_dir, i)
                         for i in self.met_subset]

        cwd = os.getcwd()
        (url, rev) = get_svn_info(cwd, self.cable_src)

        # delete local executable, copy a local copy and use taht
        local_exe = "cable"

        if os.path.isfile(local_exe):
            os.remove(local_exe)

        shutil.copy(self.cable_exe, local_exe)
        self.cable_exe = local_exe

        return (met_files, url, rev)

    def clean_up_old_files(self, site):

        out_fname = os.path.join(self.output_dir, "%s_out.nc" % (site))

        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_log_fname = os.path.join(self.log_dir, "%s_log.txt" % (site))

        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        return (out_fname, out_log_fname)

    def run_me(self, nml_fname):

        # run the model
        if self.verbose:
            os.system("./%s %s" % (self.cable_exe, nml_fname))

        else:  # No outputs to the screen, stout and stderr to dev/null
            os.system("./%s %s> /dev/null 2>&1" % (self.cable_exe, nml_fname))

    def tidy_after_run(self, met_files, url, rev):

        for fname in met_files:

            site = os.path.basename(fname).split("Fluxnet")[0]
            nml_fname = "cable_%s.nml" % (site)
            out_fname = os.path.join(self.output_dir, "%s_out.nc" % (site))

            add_attributes_to_output_file(nml_fname, out_fname, url, rev)
            shutil.move(nml_fname, os.path.join(self.namelist_dir, nml_fname))

            if self.fixed_lai is not None or self.lai_dir is not None:
                os.remove("%s_tmp.nc" % (site))


if __name__ == "__main__":

    # ------------- Change stuff ------------- #

    met_dir = "met_data"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restart_files"
    aux_dir = "CABLE-AUX/"
    namelist_dir = "namelists"
    nml_fname = "cable.nml"
    veg_fname = "def_veg_params_zr_clitt_albedo_fix.txt"
    soil_fname = "def_soil_params.txt"
    grid_fname = "gridinfo_CSIRO_1x1.nc"
    phen_fname = "modis_phenology_csiro.txt"
    cnpbio_fname = "pftlookup_csiro_v16_17tiles.csv"
    co2_conc = 37. * 10.1325
    cable_src = "CMIP6-MOSRS"
    cable_exe = "cable"
    verbose = False
    mpi = True
    num_cores = 10  # set to a number, if None it will use all cores...!

    # if empty...run all the files in the met_dir
    met_subset = []  # 'HyytialaFluxnet_met.nc']
    lai_dir = "lai_clims"
    fixed_lai = None

    # ------------------------------------------- #

    C = RunCable(met_dir, log_dir, output_dir, restart_dir, aux_dir,
                 namelist_dir, nml_fname, veg_fname, soil_fname, grid_fname,
                 phen_fname, cnpbio_fname, lai_dir, fixed_lai, co2_conc,
                 met_subset, cable_src, cable_exe, mpi, verbose)
    C.main(num_cores)

