#!/usr/bin/env bash

################################################################################
#                                                                              #
# This file is part of the TractLSM project.                                   #
#                                                                              #
# Copyright (c) 2019 Manon E. B. Sabot                                         #
#                                                                              #
# Please refer to the terms of the MIT License, which you should have received #
# along with the TractLSM.                                                     #
#                                                                              #
# This creates a text file where specific run information is set. By default,  #
# this run text file is named irun.txt and is located in the directory         #
# immediately above the TractLSM. Have a peak at this run file to see how it   #
# must be structured, or have a look below at the variables that must be       #
# specified!                                                                   #
#                                                                              #
# N.B.: if running on windows, applying the "sed -i 's/\r$//' filename"        #
#       command on the file might be necessary before it can be made an        #
#       executable.                                                            #
#                                                                              #
# __author__ = "Manon E. B. Sabot"                                             #
# __version__ = "1.0 (06.07.2018)"                                             #
# __email__ = "m.e.b.sabot@gmail.com"                                          #
#                                                                              #
################################################################################


#######################################
# Main: deals with undeclared vars,
#       writes run file
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   runfile
#######################################

main(){

# standard sanity checks and file locations
this_dir="$(cd "${BASH_SOURCE%/*}"; pwd)"

# deal with linux & windows compatibility
if [[ ${this_dir} == *\/* ]]; then
    sl="/"

elif [[ ${this_dir} == *\\* ]]; then
    sl="\\"

else
    echo "unknown dir path format"
    exit
fi

# make sure run file is located just above the TractLSM dir...
if [[ ${this_dir} == *"TractLSM"* ]]; then # if below the TractLSM
    main_dir=$(echo ${this_dir} |awk -F ${sl}"TractLSM" '{print $1}' |xargs)

    if [[ ${main_dir} == *"src"* ]]; then
        data_dir=$(echo ${main_dir} |awk -F ${sl}"src" '{print $1}' |xargs)
    fi

else
    main_dir=$(find ${this_dir} -name "TractLSM" -type d)

    if [[ ${main_dir} == "" ]]; then # the TractLSM is not present in "main dir"

        for d in $(find ${this_dir}) ; do # looking for it in dirs below

            if [[ ${d} != ${this_dir} ]]; then
                main_dir=$(find $d -name "TractLSM" -type d)

                if [[ ${main_dir} != "" ]]; then # the TractLSM is in "d"
                    main_dir=${d}
                fi
            fi

        done

        if [[ ${main_dir} == "" ]]; then # looking for it in dirs above
            if [[ ${sl} == "/" ]]; then 
                d=${this_dir%/*} # parent directory

            elif [[ ${sl} == "\\" ]]; then 
                d=${this_dir%\\*} # parent directory
            fi

            if [[ $(find ${d} -name "TractLSM" -type d) ]]; then
                 main_dir=${d}

            else
                echo "Neither above or below the TractLSM, exiting now"
                exit
            fi

        fi

    else # if the TractLSM is present in the "main dir"
        main_dir=${this_dir}
    fi

    if [[ ${main_dir} == *"src"* ]]; then
        data_dir=$(echo ${main_dir} |awk -F ${sl}"src" '{print $1}' |xargs)
    fi

fi

# default site-level files if nothing specified
if [[ ${main_dir} == *"src"* ]]; then
    runyrs="${data_dir}${sl}input${sl}fluxsites${sl}run_yrs_sites.txt"
    site_params="${data_dir}${sl}input${sl}fluxsites${sl}info_sites.csv"

else
    runyrs="${main_dir}${sl}input${sl}fluxsites${sl}run_yrs_sites.txt"
    site_params="${main_dir}${sl}input${sl}fluxsites${sl}info_sites.csv"
fi

get_undeclared_variables
write_runfile

}


#### other functions are defined in this section ###############################

#######################################
# Usage descprition for this script
# Globals:
#   None
# Arguments:
#   None
# Returns:
#   Prints usage message in shell
#######################################

usage(){ # user manual on how to use this script

cat << EOF
usage: $0 generates irun.txt, which then setups the model run(s).
None of the options are mandatory. If nothing is specified, the doy will be 1,
i.e. Jan 1st, with a run duration of 1 day; all the other default modes will
also be put in.

OPTIONS:
   -h      show this message
   -s      site name if running at site level, e.g. arg=Rocca
   -a      actual run? or not (spinup)? No arg needed
   -y      run year, e.g. arg=2003. If not set, finds match in run-yrs file
   -d      doy if running with weather generator, e.g. arg=1.
   -N      number of data days generated by the weather generator, e.g. arg=2.
   -p      name (& path) of input parameter file. Arg must be .csv or .py class.
           If not specified, but site level run, finds match in site param file 
   -l      number of time steps to shorten length of run, e.g. arg=24 is 1 day
   -m      C-H2O models to be run, e.g. arg=Control
   -r      resolution for the optimisation solver, e.g. arg=low
   -f      frequency at which the optimisation solver is called, e.g. arg=3
           calls every 3rd step of the inputted data
   -w      solving window around the previous optimisation, e.g. arg=20 binds
           the optimisation to be +/- 20% of the previous optimisation
   -t      run name tag, e.g. arg=test. The default tag is based on above -m arg
   -P      photosynthetis model, either Farquhar (default), or arg=Collatz
   -F      figure outputs using standard built-in plots? No arg needed
   -R      out repo/project name, for organisational purposes, e.g. arg=ProjectA
   -o      specific name of the run file created
EOF

}


#######################################
# Retrieve undeclared variables
# Globals:
#   None
# Arguments:
#   None
# Returns:
#   lots
#######################################

get_undeclared_variables(){

if [[ -z ${sitename+x} ]]; then # default to weather generator
    sitename="None"
    year2run="None"
    actual="yes"
    plot="yes"
fi

if [[ -z ${year2run+x} ]] && [[ $sitename != "None" ]]; then # default to 1st dr
    printf "No year supplied. Picking first year in run_yrs_sites.txt\n"

    if [[ -z ${actual+x} ]]; then # default to previous year
        actual="None" # spinup runs the previous year to that selected
    fi

    year2run=$(grep -h ${sitename} ${runyrs} |awk -F: '{print $2}' \
               |awk -F, '{print $1}' |xargs) # retrieves 1st matching year
fi

if [[ -z ${actual+x} ]]; then # default to previous year
    actual="None" # spinup runs the previous year to that selected
fi

if [[ -z ${doy+x} ]]; then # default to Jan 1st
    doy="1."
fi

if [[ -z ${Ndays+x} ]]; then # default to 1 day run
    Ndays="1."
fi

if [[ -z ${paramfile+x} ]]; then
    if [[ ${sitename} != "None" ]]; then
        paramfile=${site_params}

    else # default is the default_params class in TractLSM/Utils
        paramfile="default"
    fi
fi

if [[ -z ${Nsteps+x} ]]; then # default to full length of timeseries
    Nsteps="default"
fi

if [[ -z ${models+x} ]]; then # default is Control, ProfitMax
    models="default" 
fi

if [[ $models == "all" ]]; then
    models="Control, ProfitMax12, Wolf16"
fi

if [[ -z ${resolution+x} ]]; then # default to low resolution
    resolution="default" 
fi

if [[ -z ${fopt+x} ]]; then # default to optimising at every time step
    fopt="default" 
fi

if [[ -z ${window+x} ]]; then # default to full vulnerability curve
    window="default" 
fi

if [[ -z ${tag+x} ]]; then # default tag based on above chosen models
    tag="default" 
fi

if [[ -z ${photomod+x} ]]; then # default is the Farquhar model
    photomod="default"
fi

if [[ -z ${plot+x} ]]; then
    plot="None"
fi

if [[ -z ${project+x} ]]; then
    project="None"
fi

if [[ -z ${runfile+x} ]]; then # default name of the file created by this script
    runfile=${main_dir}${sl}"irun.txt" 
fi

}


#######################################
# Write run file
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   runfile
#######################################

write_runfile(){

touch ${runfile}

cat > ${runfile} << EOF
### The TractLSM can be used in several configurations, determined by a number
### of switches supplied to run the model. Therefore, this text file is a setup
### file equivalent to the namelist in LSMs.

# Throughout the file, 'None' or 'default' entries are the equivalent of
# commenting out the lines.
# When giving a file name, the absolute path is needed.


## 1. Site level runs

# The site can be set to match a meteorological forcings file for a specific
# site, e.g. 'Rocca'. Alternatively, 'None' calls a weather generator which uses
# idealised meteorological forcings to run the model.
# At a given site, the run year must also be specifed. If using the weather
# generator, the year must be set to 'None'.

site = $sitename
year = $year2run


# When running the model at site level, the actual switch can be set to 'no' to
# run a spinup (year -1). This spinup will be used to initialise the soil
# moisture state of the actual site and year run.
# If the switch is set to 'yes', then the given year with either the soil
# moisture initialisation from the spinup run (when it exists), or the
# parameterised saturation point, is run.

actual = $actual


## 2. Specifications for the weather generator

# The weather generator should preferably be called upon when running an
# idealised experiment. Yet, a start day (e.g. doy = 1, i.e. Jan 1st) and a
# duration for the run (e.g. Nday = 1. runs for the start day only) have to be
# chosen. If the model is being run at site level, it does not matter what the
# doy and Ndays switches equate to.

doy = $doy
Ndays = $Ndays


## 3. Using a specific set of parameters

# Specific parameters should either be kept in a csv file (e.g. containing
# sites' info) or in a class within a python file (setup like the default
# parameters in TractLSM/Utils/).
# When importing the parameters from a python class, the import can differ
# depending on your machine. So you might have to modify the retrieve_class
# function in TractLSM/Utils/general_utils.py, which is not recommended.

params = $paramfile


## 4. Choosing the model configuration

# If you do not want to run the model for the full length of the timeseries
# (e.g. only running for the first 10 hours for testing purposes), you must
# specify the number of timesteps if must be run for. If no timestep number is
# specified, the full length of the timeseries is run.

Nsteps = $Nsteps


# The models switch allows you to pick the C-H2O flux solver(s). If the switch
# is left unspecified (i.e. commented out or set to 'None'), the default
# configurations are run.
# These default choices are: 'Control, ProfitMax'.
# All possible choices are: 'Control, ProfitMax12, Wolf16'.
# Any single one of those configurations can be run alone.

models = $models


# The Farquhar photosynthesis model is the default recommended photosynthesis
# model. The Collatz model can be run instead by chosing 'Collatz'.

photo = $photomod


# The matricial optimisation solver runs at a relatively 'low' resolution, by
# default. When running the model on short timescales, it is worth exploring the
# effects of running at 'med' or 'high' resolutions.

resolution = $resolution


# By default, the optimisation solver is called for every timestep provided in
# the forcing data. However, it can also be called to run for a given frequency.
# For example, a frequency of 3 will only call the optimisation solver every 3rd
# timestep in the forcing data.

fopt = $fopt


# By default, the optimisation solver is allowed to search for the point of
# maximum profit on the full hydraulic stream. To speed up the optimisation,
# restrictions can be applied to limit the portion of the stream on which to
# search for maximum profit. For example, a window of '20' binds the search in
# +/- 20% range relative to the previous maximum profit point.
# N.B.: one must be careful not to introduce solving artefacts in combining
# this option with the fopt option

window = $window


# The tag refers to the unique name associated with the model configuration of
# the run. By default, this is based on the models that have been called upon,
# such that the tag 'CP' refers to the Control + ProfitMax.
# However, any tag can be specified by the user.

tag = $tag


## 5. Generating default figures

# If plot is set to 'yes', several figures are generated, depending on the
# configuration of the run. These can be: a comparison of the different
# optimisation configurations, a comparison of the different optimisation
# configurations with the Control, or a comparison of the Farquhar and Collatz
# photosynthesis models.

plot = $plot


## 6. Assigning the model run to a specific project

# When no project name is specified, all input and output files are stored
# directly in the input & output repositories. Otherwise, input and output files
# are stored in folders named after the project, e.g. input/var_kmax/ or
# output/var_kmax.

project = $project
EOF

sed -i 's/, */, /g' ${runfile} # make sure there is a space before each comma

}


#### execute functions #########################################################

# filter options that exist from those that don't
while getopts "hs:ay:d:N:p:l:m:r:f:w:t:P:FR:o:" OPTION; do

    case ${OPTION} in

        h)
            usage
            exit 1
            ;;

        s)
            sitename=${OPTARG}
            ;;

        a)
	        actual="yes"
	        ;;

        y)
            year2run=${OPTARG}
            ;;

        d)
            doy=${OPTARG}
            ;;

        N)
            Ndays=${OPTARG}
            ;;

        p)
            paramfile=${OPTARG}
            ;;

        l)
            Nsteps=${OPTARG}
            ;;

        m)
            models=${OPTARG}
            ;;

        r)
            resolution=${OPTARG}
            ;;

        f)
            fopt=${OPTARG}
            ;;

        w)
            window=${OPTARG}
            ;;

        t)
            tag=${OPTARG}
            ;;

        P)
            photomod=${OPTARG}
            ;;

        F)
            plot="yes"
            ;;

        R)
            project=${OPTARG}
            ;;

        o)
            runfile=${OPTARG}
            ;;

        ?)
            usage
            exit
            ;;

    esac

done

# execute main
main

