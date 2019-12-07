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
# This script can run the ProfitMax with different resolutions for the         #
# hydraulic stream, different frequencies of optimisation, and different       #
# windows of optimisation. A baseproject doesn't have to be specified, but if  #
# the forcing files exist in a different project, it's worth specifying as     #
# that will avoid running the spinup, and automatically copy the forcings to   #
# the new project.                                                             #
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
base_project="var_kmax_sample" # project name required
project="var_accuracy"  # project name required, it cannot be just input/output
#sns="Soroek7,ElSaler1k20" # sites to run (default picks best kmax calibs)
#sitenames=${sns}
resolutions="low,med" # resolution of the solve on the hydraulic stream
frequencies="default,3,6,9,24" # every 1.5, 3, and 4.5 hours
windows="default,50,200,500" # restricted to x% around morning stream

# NB1: to avoid ALL combinations of frequencies and windows from being run,
# "default" must be specified first in those lists.

# NB2: the high resolution is really high and very long to run, so using it is
# not recommended. "low" is recommended for testing purposes; "med" seems
# sufficient for increased accuracy on short timescales.



#######################################
# Main: generates runfiles, runs model
#       at multiple sites with
#       variations on the solving
#       resolution, frequency, etc.
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   input, output, & log files
#######################################

main(){

# default site met input rep
met="${data_dir}${sl}input${sl}fluxsites${sl}met${sl}"
suffix="Fluxnet_met.nc" # suffix common to all the met files

# run file generator
mkrunfile="${this_dir}${sl}preparation_scripts${sl}create_setup_file.sh"

# model run launcher
runmodel="${this_dir}${sl}ExecTractLSM"

get_undeclared_variables # retrieve variables if undeclared
get_combinations # get all combinations of site-year to run
initialise_logfile # write the headers in the log file

i=0 # runfile number tracker
stored_combis=() # combis storer, to run actual after spinup
child_process_ids="" # track job ids to check for successful completion
child_subprocess_ids="" # track subjob ids to stall until completed

for combis in "${combinations[@]}"; do # loop over combinations

    sitename=$(echo ${combis} |awk -F: '{print $1}')
    year2run=$(echo ${combis} |awk -F: '{print $2}')
    rez=$(echo ${combis} |awk -F: '{print $3}')
    freq=$(echo ${combis} |awk -F: '{print $4}')
    win=$(echo ${combis} |awk -F: '{print $5}')

    if [[ ${project} != "None" ]]; then
        if [[ ! -d  "${data_dir}${sl}input${sl}projects${sl}${project}" ]]; then
            mkdir "${data_dir}${sl}input${sl}projects${sl}${project}"
        fi

        if [[ ! -d  "${data_dir}${sl}output${sl}projects${sl}${project}" ]]
        then
            mkdir "${data_dir}${sl}output${sl}projects${sl}${project}"
        fi

        if [[ ! -z "${base_project}" ]]; then # cp forcings from baseproject
            inputfile1="${data_dir}${sl}input${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}input${sl}projects${sl}${project}"
            inputfile1+="${sl}${sitename}_met_and_plant_data_spinup_"
            inputfile2+="${sl}${sitename}_met_and_plant_data_spinup_"
            inputfile1+="${year2run}.csv"
            inputfile2+="${year2run}.csv"
            cp "${inputfile1}" "${inputfile2}"

            inputfile1="${data_dir}${sl}output${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}output${sl}projects${sl}${project}${sl}"
            inputfile1+="${sl}${sitename}_"
            cp "${inputfile1}"*"spinup_${year2run}.csv" "${inputfile2}"

            inputfile1="${data_dir}${sl}input${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}input${sl}projects${sl}${project}"
            inputfile1+="${sl}${sitename}_met_and_plant_data_actual_"
            inputfile2+="${sl}${sitename}_met_and_plant_data_actual_"
            inputfile1+="${year2run}.csv"
            inputfile2+="${year2run}.csv"
            cp "${inputfile1}" "${inputfile2}"

            inputfile1="${data_dir}${sl}output${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}output${sl}projects${sl}${project}${sl}"
            inputfile1+="${sl}${sitename}_"
            cp "${inputfile1}"*"actual_${year2run}.csv" "${inputfile2}"
            wait
        fi

    else
        if [[ ! -z "${base_project}" ]]; then # cp forcings from baseproject
            inputfile1="${data_dir}${sl}input${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}input"
            inputfile1+="${sl}${sitename}_met_and_plant_data_spinup_"
            inputfile2+="${sl}${sitename}_met_and_plant_data_spinup_"
            inputfile1+="${year2run}.csv"
            inputfile2+="${year2run}.csv"
            cp "${inputfile1}" "${inputfile2}"

            inputfile1="${data_dir}${sl}output${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}output"
            inputfile1+="${sl}${sitename}_"
            cp "${inputfile1}"*"spinup_${year2run}.csv" "${inputfile2}"

            inputfile1="${data_dir}${sl}input${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}input${sl}"
            inputfile1+="${sl}${sitename}_met_and_plant_data_actual_"
            inputfile2+="${sl}${sitename}_met_and_plant_data_actual_"
            inputfile1+="${year2run}.csv"
            inputfile2+="${year2run}.csv"
            cp "${inputfile1}" "${inputfile2}"

            inputfile1="${data_dir}${sl}output${sl}projects${sl}${base_project}"
            inputfile2="${data_dir}${sl}output${sl}"
            inputfile1+="${sl}${sitename}_"
            cp "${inputfile1}"*"actual_${year2run}.csv" "${inputfile2}"
            wait
        fi
    fi

    runrun "spinup" &
    child_subprocess_ids+="$! " # store subjob ids
    stored_combis+=("${combis}") # store combis
    let i+=1 # update runfile number tracker

    # too many combi for allowed # of CPUs or last combi / element
    if [[ (( "${i}" == "${NCPUs}" )) || \
          (( "${combis}" == "${combinations[${#combinations[@]}-1]}" )) ]]
    then

        # wait for subjobs to finish, then continue
        for child_subprocess_id in ${child_subprocess_ids}; do

            { wait ${child_subprocess_id} ; }

        done

        child_process_ids+="${child_subprocess_ids}" # subprocesses
        child_subprocess_ids="" # reset subjob id tracker
        i=0 # reset runfile number tracker

        for recombis in "${stored_combis[@]}"; do

            sitename=$(echo ${recombis} |awk -F: '{print $1}')
            year2run=$(echo ${recombis} |awk -F: '{print $2}')
            rez=$(echo ${recombis} |awk -F: '{print $3}')
            freq=$(echo ${recombis} |awk -F: '{print $4}')
            win=$(echo ${recombis} |awk -F: '{print $5}')
            runrun "actual" &
            child_subprocess_ids+="$! " # store subjob ids
            let i+=1 # update runfile number tracker

        done

        # wait for subjobs to finish, then continue
        for child_subprocess_id in ${child_subprocess_ids}; do

            { wait ${child_subprocess_id} ; }

        done

        child_process_ids+="${child_subprocess_ids}" # subprocesses
        child_subprocess_ids="" # reset subjob id tracker
        stored_combis=() # reset combis storer
        i=0 # reset runfile number tracker
    fi

done

get_runstatus # check whether the runs were successful

# IOError (e.g. no such file or directory) in log, no tidy
if ! grep -q "IOError: " ${log_file}; then
    tidy # remove unnecessary run files created
fi

}


#### other functions are defined in this section ###############################

#######################################
# Default to all sites
# Globals:
#   met, suffix
# Arguments:
#   None
# Returns:
#   array containing all site names
#######################################

all_sites(){

sites=()

for filename in ${met}*.nc; do

    sitename=$(echo $(basename -- "${filename}") |awk -F${suffix} '{print $1}')
    sites+=("${sitename}")

done

}


#######################################
# Default to all years
# Globals:
#   sites, runyrs
# Arguments:
#   None
# Returns:
#   array containing site run years
#######################################

all_years(){

years2run=()

for sitename in "${sites[@]}"; do # retrieve years

    yearsofsite=$(grep -h ${sitename} ${runyrs} |awk -F: '{print $2}' \
                  |sed 's/\, /-/g' |xargs)

    if [[ -z "${yearsofsite}" ]]; then # if sitename has 1 character extension
        yearsofsite=$(grep -h "${sitename%?}" ${runyrs} |awk -F: '{print $2}' \
                      |sed 's/\, /-/g' |xargs)
    fi

    if [[ -z "${yearsofsite}" ]]; then # if sitename has 2 character extension
        yearsofsite=$(grep -h "${sitename%??}" ${runyrs} |awk -F: '{print $2}' \
                      |sed 's/\, /-/g' |xargs)
    fi

    if [[ -z "${yearsofsite}" ]]; then # if sitename has 3 character extension
        yearsofsite=$(grep -h "${sitename%???}" ${runyrs} \
                      |awk -F: '{print $2}' |sed 's/\, /-/g' |xargs)
    fi

    years2run+=("${yearsofsite}")

done

}


#######################################
# Retrieve undeclared variables
# Globals:
#   met, suffix, sites, runyrs
# Arguments:
#   None
# Returns:
#   sites, years2run,
#   (models, project)
#######################################

get_undeclared_variables(){

# NCPUs for parallelisation
if [[ -z ${NCPUs+x} ]]; then 
    NCPUs=$(grep ^cpu\\scores /proc/cpuinfo |uniq |awk '{print $4}')
    NCPUs=$((${NCPUs} / 2)) # only take half the CPUs
fi

if [[ -z ${models+x} ]]; then # default is Control, ProfitMax
    models="default" 
fi

# is the project specified?
if [[ -z ${project+x} ]]; then
    project="None"
    runyrs="${data_dir}${sl}input${sl}fluxsites${sl}run_yrs_sites.txt"

else
    runyrs="${data_dir}${sl}input${sl}projects${sl}${project}${sl}run_yrs_sites.txt"

    if [[ ! -f ${runyrs} ]]; then
        runyrs="${data_dir}${sl}input${sl}fluxsites${sl}run_yrs_sites.txt"
    fi
fi

# default to all sites in input/fluxsites/met
if [[ -z ${sites+x} ]]; then 
      all_sites
fi

# default to all years in run_yrs_sites.txt
if [[ ( ( -n ${years2run+x} ) && ( ${#years2run[@]} != ${#sites[@]} && \
          ${#years2run[@]} != 1 ) ) || ( -z ${years2run+x} ) ]]; then
    all_years
fi

# is the photosynthesis model specified?
if [[ -z ${photomod+x} ]]; then
    photomod="default"
fi

}


#######################################
# Get all combinations of site-year
# Globals:
#   sites, years2run
# Arguments:
#   None
# Returns:
#   combinations of site-year
#######################################

get_combinations(){

i=0 # initialise yearsofsite tracker
combinations=()

for ss in "${sites[@]}"; do # loop over sites

    for yearsofsite in "${years2run[i]}"; do # loop over years of all sites

        # make string of years into array
        yearsofsite=$(echo ${yearsofsite} |sed 's,\-, ,g')
        yearsofsite=(${yearsofsite// / })

        for yy in "${yearsofsite[@]}"; do # loop over years of this site

            for rez in $(echo ${resolutions} |sed "s/,/ /g"); do

                for freq in $(echo ${frequencies} |sed "s/,/ /g"); do

                    for win in $(echo ${windows} |sed "s/,/ /g"); do

                        if [[ (( ${rez} == "low" )) && \
                              (( ${freq} == "default" )) && \
                              (( ${win} == "default" )) ]]; then
                            :

                        elif [[ (( ${rez} == "med" )) && \
                                (( ${freq} == "default" )) && \
                                (( ${win} == "default" )) ]]; then

                            combinations+=("${ss}:${yy}:${rez}:${freq}:${win}")

                        elif [[ (( ${rez} == "high" )) && \
                                (( ${freq} == "default" )) && \
                                (( ${win} == "default" )) ]]; then

                            combinations+=("${ss}:${yy}:${rez}:${freq}:${win}")

                        elif [[ (( ${rez} == "low" )) && \
                                (( (( ${freq} == "default" )) || \
                                   (( ${win} == "default" )) )) ]]; then

                            combinations+=("${ss}:${yy}:${rez}:${freq}:${win}")
                        fi

                    done

                done

            done

        done

    done

    let i+=1 # update yearsofsite tracker

done

}


#######################################
# Write headers in log file
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   runfile
#######################################

initialise_logfile(){

if [[ $(echo $(ps -o stat= -p $PPID)) == *"s"* ]]; then
    if [[ $(echo $(ps -o stat= -p $$)) == *"+"* ]]; then
        launloc=$(echo ${BASH_SOURCE[0]} |tail -c +3)

    else
        launloc=$(ps -f -p $PPID |tail -1 |awk '{print $NF}'|tail -c +3)
    fi

else
    launloc=$(ps -f -p $PPID |tail -1 |awk '{print $NF}'|tail -c +3)
fi

cat > ${log_file} << EOF
### The log file was started on $(date +'%d/%m/%Y %H:%M:%S')
###
### The project is: ${project}
###
### The sites ran are: ${sitenames[@]}
### The years ran are: ${years2run[@]}
###
### The models ran are: ${models}
### The photosynthesis models is: ${photomod}
### The hydraulic stream's resolutions are: ${resolutions}
### The solving frequencies are: ${frequencies}
### The solving windows are: ${windows}
###
### These model runs were launched from ${launloc}
###

EOF

sed -i 's/, */, /g' ${log_file} # make sure there is a space before each comma

}


#######################################
# Run spinup, or actual year
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   input files & run outputs
#######################################

runrun(){

runfile="${this_dir}${sl}tmp${sl}runme${i}.txt"

if [[ ${1} == "spinup" ]]; then # make irunN.txt for spinup & run model
    ${mkrunfile} -s "${sitename}" -y "${year2run}" -m Control -P "${photomod}" \
                 -R "${project}" -o "${runfile}" >> ${log_file} 2>&1

elif [[ ${1} == "actual" ]]; then # make irunN.txt for actual & run model
    ${mkrunfile} -s "${sitename}" -a -y "${year2run}" -m "${models}" \
                 -r "${rez}" -f "${freq}" -w "${win}" -P "${photomod}" \
                 -R "${project}" -o "${runfile}" >> ${log_file} 2>&1
fi

if [[ ${1} == "spinup" ]]; then
    echo "Running spinup" ${sitename} ${year2run}
    echo "The ${sitename} ${year2run} spinup run" >> ${log_file}

elif [[ ${1} == "actual" ]]; then
    echo "Running actual" ${sitename} ${year2run} ${rez} ${freq} ${win}
    message="The ${sitename} ${year2run} ${photomod} ${rez} ${freq} ${win} "
    message+="actual run"
    echo "${message}" >> ${log_file}   
fi

${runmodel} ${runfile} >> ${log_file} 2>&1

}


#######################################
# Check success of subprocesses
# Globals:
#   child_process_ids, combinations
# Arguments:
#   None
# Returns:
#   success or failure message
#######################################

get_runstatus(){

i=3 # text file info on the specific runs starts at line 3

for child_process_id in ${child_process_ids}; do

    # pass each process id into the wait command to retrieve its exit code
    { wait ${child_process_id} ; }
    rc=$? # store exist code in rc

    # check that specific line is not an error message, if it is, move to next
    while true; do

        if [[ ( $(sed -n "${i}p" < ${log_file}) == *"The"* ) && \
              ( $(sed -n "${i}p" < ${log_file}) == *"run"* ) ]]; then
            break

        else
            let i+=1 # go to next line
        fi

    done

    # inspect each process exit codes
    if [ $rc -ne 0 ]; then
        sed "${i}s/$/ failed with the exit code $rc/" -i ${log_file}

    else
        sed "${i}s/$/ was successful/" -i ${log_file}
    fi

    let i+=1 # update line number

done

}


#######################################
# Tidy unnecessary created files
# Globals:
#   None
# Arguments:
#   None
# Returns:
#   None
#######################################

tidy(){

to_rm=$(find ${this_dir}${sl}tmp${sl}"runme"*[0-9]*".txt" \
        |grep -v "${this_dir}${sl}tmp${sl}runme.txt" |grep -v ${log_file} \
        |xargs)
rm ${to_rm}

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

# find the best calibs if the site names have not been specified
if [[ -z ${sitenames+x} ]]; then
    cals="${data_dir}${sl}output${sl}projects${sl}"
    cals+="${base_project/sample/best_calib}${sl}best.txt"
    sitenames=""

    while IFS= read -r line; do

        sitenames+=$(echo "${line//: },")

    done < ${cals}

    sitenames=$(echo "${sitenames::-1}")
fi

# make site names into list
sites=()

for i in $(echo ${sitenames} |sed "s/,/ /g"); do

    sites+=("${i}")

done

# log file
log_file="${this_dir}${sl}tmp${sl}log.o$(date +%s)"

# execute main & output run time in log
{ time main ; } >> ${log_file} 2>&1
wait
