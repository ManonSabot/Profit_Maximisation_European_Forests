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
# This script applies user defined operations on the forcing data and then     #
# proceeds to rerun the model with those changes. It can be useful in order to #
# perform sensitivity analyses or to calibrate a variable. Here, it is         #
# hardwired to change the soil moisture stress function in the control, #
# but this can easily be altered.                                              #
#                                                                              #
# N.B.: if running in windows, applying the "sed -i 's/\r$//' filename"        #
#       command on the file might be necessary before it can be made an        #
#       executable.                                                            #
#                                                                              #
# __author__ = "Manon Sabot"                                                   #
# __version__ = "1.0 (23.10.2019)"                                             #
# __email__ = "m.e.b.sabot@gmail.com"                                          #
#                                                                              #
################################################################################


#### USER'S CHOICES ####
NCPUs=4
models="Control"
base_project="calib_fw" # input in based project to modify
runyrs="run_yrs_sites.txt"
sns="Soroe,ElSaler1,Hyytiala,Hesse,Puechabon,Parco,Rocca1,Rocca2,Loobos,Espirra"
sites=${sns}


#######################################
# Main: runs model at multiple sites
#       for variations on a given
#       model param
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   input, output
#######################################

main(){

project=${base_project}

# if no base directory or files in there, first run the default config
if [[ ! -d ${data_dir}${sl}input${sl}projects${sl}${project} ]]
    then
        mkdir ${data_dir}${sl}input${sl}projects${sl}${project}
fi

if [[ ! -d ${data_dir}${sl}output${sl}projects${sl}${project} ]]
    then
        mkdir ${data_dir}${sl}output${sl}projects${sl}${project}
fi

if [[ ! "$(ls ${data_dir}${sl}output${sl}projects${sl}${project}${sl}*.csv)" ]]
    then
        echo "initialising the reference input/output..."
        ${this_dir}${sl}run_multiple_sites.sh -c ${NCPUs} -m $models \
                                              -R $project -s $sites
        wait # make sure the job is finished before moving on
        echo "initialisation performed, now moving on"
fi

# now optimise the parameters
if [[ ! -z ${ssl} ]]; then 
    sites=${ssl} # temporarily change the site list
fi

sample

if [[ ! -z ${ssl} ]]; then
    sites=${sns} # change the site list back
fi

}


#### other functions are defined in this section ###############################

#######################################
# Input files for the various fws
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   input files
#######################################

var_fws(){

# loop over the site names
for i in $(echo ${sites} |sed "s/,/ /g"); do

    # create the new actual sites input files with updated fws
    wf=${data_dir}${sl}input${sl}projects${sl}${project}${sl}
    wf+=${i}*"actual"*".csv"
    filenames=$(find ${wf})

    for j in ${filenames[@]}; do

        ${this_dir}${sl}preparation_scripts${sl}parameter_sampling.py "${j}" \
                                                                      -s "fw"

    done

done

}


#######################################
# Sample the model parameters &
# run the model for all sample values
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   data files
#######################################

sample(){

# run all configs with adjustment
if [[ ! -d ${data_dir}${sl}input${sl}projects${sl}${base_project} ]]; then
    echo "The project the adjustment is based on isn't present, fix that!"
    exit 1
fi

var_fws

# retrieve all new site names
all_sites=${data_dir}${sl}input${sl}projects${sl}${project}${sl}*"actual"*".csv"
all_sites=($all_sites)

for i in "${!all_sites[@]}"; do

    all_sites[${i}]="${all_sites[${i}]##*/}"
    all_sites[${i}]=$(echo ${all_sites[${i}]} |sed 's/_.*//g')
    all_sites[${i}]="${all_sites[${i}]%%_met_and_plant_data_actual}"

done

# make sure the site list doesn't replicate names (if multiple actual years)
all_sites=($(echo "${all_sites[@]}" |tr ' ' '\n' |sort -u |tr '\n' ' '))

# convert array to list
all_sites=$( IFS=$','; echo "${all_sites[*]}" )

${this_dir}${sl}run_multiple_sites.sh -c ${NCPUs} -m $models -R $project \
                                      -s $all_sites

}


#######################################
# Assess the performance of the
# various model calibrations
# Globals:
#   lots
# Arguments:
#   adjustment type (optional)
# Returns:
#   data files
#######################################

performance(){

# clean up existing best files
rm ${data_dir}${sl}output${sl}projects${sl}${base_project}${sl}best.txt

cal=${data_dir}${sl}src${sl}postprocessing_scripts${sl}parameter_calibration.py

for i in $(echo ${sites} |sed "s/,/ /g"); do

    python ${cal} ${i} ${base_project} -p "fw"

done

}

#### execute functions #########################################################

# standard sanity checks and file locations
this_dir="$(cd "${BASH_SOURCE%/*}"; pwd)"

# deal with linux & windows compatibility
if [[ ${this_dir} == *\/* ]]; then
    sl="/"

elif [[ ${this_dir} == *\\* ]]; then
    sl="\\"

else
    echo "Unknown dir path format, now exiting."
    exit 1
fi

# make sure this is located above the TractLSM dir...
if [[ ${this_dir} == *"TractLSM"* ]]; then
    this_dir=$(echo ${this_dir} |awk -F ${sl}"TractLSM" '{print $1}' |xargs)
fi

if [[ ${this_dir} == *"src"* ]]; then
    data_dir=$(echo ${this_dir} |awk -F ${sl}"src" '{print $1}' |xargs)

else
    data_dir=${this_dir}
fi

# make sure a reasonable number of cores has been requested
unicore=$(grep ^cpu\\scores /proc/cpuinfo |uniq |awk '{print $4}')

# check whether there is REAL hyperthreading &/or too many requested cores
if [[ (( (( ${NCPUs} -gt ${unicore} )) && (( ${unicore} -lt \
            $(grep -E "cpu cores|siblings|physical id" /proc/cpuinfo \
              |xargs -n 11 echo |sort |uniq |awk '/^physical id/{print $7}' \
              |tail -1) )) )) || (( ${NCPUs} -ge $(nproc --all) )) ]]; then

    if [[ ${NCPUs} -ge $(nproc --all) ]]; then
        let NCPUs-=1

    else
        NCPUs=${unicore}
    fi
fi

# execute main
main
wait # make sure the job is finished before moving on

# establish the best performing runs
performance
wait # make sure the job is finished before moving on

