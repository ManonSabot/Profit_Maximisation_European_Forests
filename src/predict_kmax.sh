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
# This script applies variations to the values of kmax and runs the model      #
# accordingly. Precisely, it first establishes a baseline for the long term    #
# optimisation of kmax at 25 deg C and 1 kPa (in base_project).                #
# It can proceed to optimise kmax for an average climate (in _adjust_average)  #
# and for an extreme climate (in _adjust_extreme), while also testing a range  #
# of different long term behaviours (siteH for "high" at P12, site for         #
# "optimal", siteL for "low" right before the cost offsets the net profit. See #
# calc_site_kmax.py for further detail.                                        #
# It can also sample the values of kmax around those given in the base_project #
# (in _sample), like a calibration exercise (see calc_site_kmax.py).           #
#                                                                              #
# N.B.: if running in windows, applying the "sed -i 's/\r$//' filename"        #
#       command on the file might be necessary before it can be made an        #
#       executable.                                                            #
#                                                                              #
# __author__ = "Manon E. B. Sabot"                                             #
# __version__ = "1.0 (23.08.2018)"                                             #
# __email__ = "m.e.b.sabot@gmail.com"                                          #
#                                                                              #
################################################################################


#### USER'S CHOICES ####
NCPUs=4
models="Control,ProfitMax"
base_project="var_kmax" # project name required, it cannot be just input/output
runyrs="run_yrs_sites.txt"
sns="Soroe,ElSaler1,Hyytiala,Hesse,Puechabon,Parco,Rocca1,Rocca2,Loobos,Espirra"
sites=${sns}


#######################################
# Main: runs model at multiple sites
#       for variations on the optimal
#       value of kmax
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   input, output & perfile & summary
#######################################

main(){

# read in the user's input
options

project=${base_project}

# run 'optimal kmax' / default config
${this_dir}${sl}run_multiple_sites.sh -c ${NCPUs} -m $models -R $project \
                                      -s $sites
wait # make sure the job is finished before moving on

if [[ ${variation} == "climate" ]] || [[ ${variation} == "all" ]]; then

    if [[ ! -d ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures ]]
    then
        mkdir ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures
    fi

    # generate the sitesH & sitesL lists
    sitesH=$(echo ${sites} |sed 's/,/H,/g')'H'
    sitesL=$(echo ${sites} |sed 's/,/L,/g')'L'

    # first step is to check that the site-specific broader climate info exists
    python ${this_dir}${sl}TractLSM${sl}Utils${sl}cru_climate_lat_lon.py -VPD \
                                                                         -Tair
    python ${this_dir}${sl}TractLSM${sl}Utils${sl}cru_climate_lat_lon.py -VPD \
                                                                         -Txx
    python ${this_dir}${sl}TractLSM${sl}Utils${sl}cru_climate_lat_lon.py -LAI

    # average adjustment, all three behaviours
    adjust

    if [[ ! -d ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures ]]
    then
        mkdir ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures
    fi

    # adjustment to extremes, all three behaviours
    adjust "extreme"

    if [[ ! -d ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures ]]
    then
        mkdir ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures
    fi

fi

if [[ ${variation} == "sample" ]] || [[ ${variation} == "all" ]]; then

    if [[ ! -z ${ssl} ]]; then
        sites=${ssl} # temporarily change the site list
    fi

    sample

    if [[ ! -d ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures ]]
    then
        mkdir ${data_dir}${sl}output${sl}projects${sl}${project}${sl}figures
    fi

    if [[ ! -z ${ssl} ]]; then
        sites=${sns} # change the site list back
    fi

fi

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
usage: $0 runs the model over fluxnet sites for set years and performs a series
of operations to vary how kmax gets parameterised.
The number of CPUs to parallelise things on, the model names, reference project,
site list, and run-yr file are all set by the user INSIDE this script.
Neither one of the options is mandatory.
If none is parsed, then the variation of kmax will depend on behaviour and on
two sets of climatological values of both VPD and temperature.
If -S in parsed, then the variation around the optimal value of kmax will be
like a variation around kmax,opt.
If -A is parsed, then both the "none" option and -S parsing will occur.

OPTIONS:
   -h      show this message
   -S      variation of kmax around kmax,opt
   -A      both behavioural and climatological variations of kmax & -s option
EOF

}


#######################################
# Check whether options were declared
# Globals:
#   parsed options
# Arguments:
#   None
# Returns:
#   variation
#######################################

options(){

if [[ -z ${variation} ]]; then
    variation="climate"
fi

}


#######################################
# Input files for high & low kmax &
# for certain climatologies
# Globals:
#   lots
# Arguments:
#   adjustment type (optional)
# Returns:
#   input files
#######################################

var_kmax(){

if [[ -n ${1} ]]; then # set up new repo!

    if [[ ! -d ${data_dir}${sl}input${sl}projects${sl}${base_project}_${1} ]]
    then
        mkdir ${data_dir}${sl}input${sl}projects${sl}${base_project}_${1}
    fi

    if [[ ! -d ${data_dir}${sl}output${sl}projects${sl}${base_project}_${1} ]]
    then
        mkdir ${data_dir}${sl}output${sl}projects${sl}${base_project}_${1}
    fi

    # copy input run-year info file
    cp ${data_dir}${sl}input${sl}projects${sl}${base_project}${sl}${runyrs} \
       ${data_dir}${sl}input${sl}projects${sl}${base_project}_${1}${sl}

    for i in $(echo ${sites} |sed "s/,/ /g"); do

        # copy both basic (i.e. optimal behaviour) spinup and actual input
        wf=${data_dir}${sl}input${sl}projects${sl}${base_project}${sl}${i}_*
        filenames=$(find ${wf})

        for j in ${filenames[@]}; do

            newname=$(echo "${j}" \
                      |sed 's/'"${base_project}"'/'"${base_project}_${1}"'/')

            if [[ ! -f ${newname} ]]; then
                cp ${j} ${newname}
            fi


        done

        # copy both basic (i.e. optimal behaviour) spinup and actual input
        wf=${data_dir}${sl}output${sl}projects${sl}${base_project}${sl}

        if [[ ${1} == "sample" ]]; then
            wf+=${i}_*".csv"

        else
            wf+=${i}_*"spinup"*".csv"
        fi

        filenames=$(find ${wf})

        for j in ${filenames[@]}; do

            newname=$(echo "${j}" \
                      |sed 's/'"${base_project}"'/'"${base_project}_${1}"'/')

            if [[ ! -f ${newname} ]]; then
                cp ${j} ${newname}
            fi

        done

    done

    project="${base_project}_${1}"
fi

# loop over the site names
for i in $(echo ${sites} |sed "s/,/ /g"); do

    # create the new actual sites input files with updated kmax
    wf=${data_dir}${sl}input${sl}projects${sl}${project}${sl}${i}_*"actual"*
    filenames=$(find ${wf})

    for j in ${filenames[@]}; do

        if [[ ${1} == "sample" ]]; then
            if [[ ( ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'k0_/') ) && \
                  ( ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'k23_/') ) ]]
            then
                ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                              "${j}" -b "sample"
            fi

        elif [[ ${1} == *"average"* ]]; then
            if [[ ( ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'H_/') ) && \
                  ( ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'L_/') ) ]]
            then
                ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                      "${j}" -b "all" -VPD -Tair
                mv $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'O_/') ${j} # optimal

            else
                if [[ ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'H_/') ]]
                then
                    ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                  "${j}" -b "optimal" -VPD -Tair
                    mv $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'O_/') ${j}
                    ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                     "${j}" -b "high" -VPD -Tair
                fi

                if [[ ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'L_/') ]]
                then
                    ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                      "${j}" -b "low" -VPD -Tair
                fi
            fi

        elif [[ ${1} == *"extreme"* ]]; then
            if [[ ( ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'H_/') ) && \
                  ( ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'L_/') ) ]]
            then
                ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                       "${j}" -b "all" -VPD -Txx
                mv $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'O_/') ${j} # optimal

            else
                if [[ ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'H_/') ]]
                then
                    ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                   "${j}" -b "optimal" -VPD -Txx
                    mv $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'O_/') ${j}
                    ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                      "${j}" -b "high" -VPD -Txx
                fi

                if [[ ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'L_/') ]]
                then
                    ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                       "${j}" -b "low" -VPD -Txx
                fi
            fi

        else
            if [[ ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'H_/') ]]; then
                ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                                "${j}" -b "high"
            fi

            if [[ ! -f $(echo "${j}" |sed 's/'"${i}"'_/'"${i}"'L_/') ]]; then
                ${this_dir}${sl}preparation_scripts${sl}calc_site_kmax.py \
                                                                 "${j}" -b "low"
            fi
        fi

    done

done

}


#######################################
# Retrieve kmax for a specific climate
# & run the model for all three high,
# opt, low and behaviours
# Globals:
#   lots
# Arguments:
#   adjustment type (optional)
# Returns:
#   data files
#######################################

adjust(){

# run all configs with adjustment
if [[ ! -d ${data_dir}${sl}input${sl}projects${sl}${base_project} ]]; then
    echo "The project the adjustment is based on isn't present, fix that!"
    exit 1
fi

# retrieve all behavioural adjusted kmax
if [[ -n ${1} ]]; then
    var_kmax "adjust_${1}"

else
    var_kmax "adjust_average"
fi

${this_dir}${sl}run_multiple_sites.sh -c ${NCPUs} -m $models -R $project \
                                      -s $sites
${this_dir}${sl}run_multiple_sites.sh -c ${NCPUs} -m $models -R $project \
                                      -s $sitesH
${this_dir}${sl}run_multiple_sites.sh -c ${NCPUs} -m $models -R $project \
                                      -s $sitesL

}


#######################################
# Sample kmax around a reasonable kmax
# value at Tair = 25 deg C, VPD = 1 kPa
# & run the model for all sample values
# Globals:
#   lots
# Arguments:
#   adjustment type (optional)
# Returns:
#   data files
#######################################

sample(){

# run all configs with adjustment
if [[ ! -d ${data_dir}${sl}input${sl}projects${sl}${base_project} ]]; then
    echo "The project the adjustment is based on isn't present, fix that!"
    exit 1
fi

var_kmax "sample"

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
# respective model configurations
# Globals:
#   lots
# Arguments:
#   adjustment type (optional)
# Returns:
#   data files
#######################################

performance(){

# first clean up existing files
rm ${data_dir}${sl}output${sl}projects${sl}${base_project}${sl}perf_scores.csv
rm ${data_dir}${sl}output${sl}projects${sl}${base_project}_best_climate${sl}*
rm ${data_dir}${sl}output${sl}projects${sl}${base_project}_best_calib${sl}*
rmdir ${data_dir}${sl}output${sl}projects${sl}${base_project}_best_climate
rmdir ${data_dir}${sl}output${sl}projects${sl}${base_project}_best_calib

# loop over the site names
for i in $(echo ${sites} |sed "s/,/ /g"); do

    wf=${data_dir}${sl}input${sl}projects${sl}${base_project}${sl}
    wf+=${i}_*"actual"*".csv"
    filenames=$(find ${wf})

    ${this_dir}${sl}postprocessing_scripts${sl}performance_scores.py \
        $(echo ${filenames} |awk '{print $1}')

done

}


#### execute functions #########################################################

# filter options that exist from those that don't
while getopts "hSA" OPTION; do

    case ${OPTION} in

        h)
            usage
            exit 1
            ;;

        S)
            variation="sample"
            ;;

        A)
            variation="all"
            ;;

        ?)
            usage
            exit 1
            ;;

    esac

done


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

# establish the performance scores & best performing runs
performance
wait # make sure the job is finished before moving on

# summarise the input and output
if [[ ${variation} == "sample" ]]; then
    ${this_dir}${sl}postprocessing_scripts${sl}summary_tables.py \
        ${base_project}_sample

elif [[ ${variation} == "all" ]]; then
    ${this_dir}${sl}postprocessing_scripts${sl}summary_tables.py \
        ${base_project}_adjust_average ${base_project}_adjust_extreme \
        ${base_project}_sample

else
    ${this_dir}${sl}postprocessing_scripts${sl}summary_tables.py \
        ${base_project}_adjust_average ${base_project}_adjust_extreme
fi
