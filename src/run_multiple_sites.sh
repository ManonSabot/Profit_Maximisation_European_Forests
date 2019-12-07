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
# This script either takes in site names and years for the model to be run     #
# over fluxnet site data, or it looks for the available site data before       #
# running the TractLSM over all available sites and years (see the             #
# run_yrs_sites.txt file in input/fluxsites for those), first for the spinup   #
# spinup year run, then the actual year run.                                   #
#                                                                              #
# N.B.: if running in windows, applying the "sed -i 's/\r$//' filename"        #
#       command on the file might be necessary before it can be made an        #
#       executable.                                                            #
#                                                                              #
# __author__ = "Manon E. B. Sabot"                                             #
# __version__ = "1.0 (12.07.2018)"                                             #
# __email__ = "m.e.b.sabot@gmail.com"                                          #
#                                                                              #
################################################################################


#######################################
# Main: generates runfiles, runs model
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

    sitename=$(echo ${combis}|awk -F: '{print $1}')
    year2run=$(echo ${combis}|awk -F: '{print $2}')
    runrun "spinup" &
    child_subprocess_ids+="$! " # store subjob ids
    stored_combis+=("${combis}") # store combis
    let i+=1 # update runfile number tracker

    # too many combi for allowed # of CPUs or last combi / element
    if [[ (( "${i}" == "${NCPUs}" )) || \
          (( "${combis}" == "${combinations[${#combinations[@]}-1]}" )) ]]; then
        # wait for subjobs to finish, then continue

        for child_subprocess_id in ${child_subprocess_ids}; do

            { wait ${child_subprocess_id} ; }

        done

        child_process_ids+="${child_subprocess_ids}" # store the subprocesses
        child_subprocess_ids="" # reset subjob id tracker
        i=0 # reset runfile number tracker

        for recombis in "${stored_combis[@]}"; do

            sitename=$(echo ${recombis}|awk -F: '{print $1}')
            year2run=$(echo ${recombis}|awk -F: '{print $2}')
            runrun "actual" &
            child_subprocess_ids+="$! " # store subjob ids
            let i+=1 # update runfile number tracker

        done

        # wait for subjobs to finish, then continue
        for child_subprocess_id in ${child_subprocess_ids}; do

            { wait ${child_subprocess_id} ; }

        done

        child_process_ids+="${child_subprocess_ids}" # store the subprocesses
        child_subprocess_ids="" # reset subjob id tracker
        stored_combis=() # reset combis storer
        i=0 # reset runfile number tracker
    fi

done

get_runstatus # check whether the runs were successful

# if IOError (e.g. no such file or directory) in log, do not tidy repo
if ! grep -q "IOError: " ${log_file}; then
    tidy # remove unnecessary run files created for the parallelisation
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
usage: $0 runs the model over fluxnet sites for set years.
Neither one of the two options (-s & -y) are mandatory.
A single year can be specified, even if no sites are. That will run the TractLSM
for a single year across all the sites where that year of input is available.  
If no option is specified at all, this script will look for the available met
site data before running the TractLSM over all available sites and years for
each site, as specified in input/fluxsites/run_yrs_sites.txt, first for the
spinup year run, then the actual year run.

In the example below, the co-occuring site-year are:
Rocca 2003
Espirra 2004 2005 2006

OPTIONS:
   -h      show this message
   -c      number of CPUs for parallelisation
   -m      C-H2O models to be run, e.g. arg=Control
   -P      photosynthetis model, either Farquhar (default), or arg=Collatz
   -R      out repo/project name, for organisational purposes, e.g. arg=ProjectA
   -s      site name(s) if running at site level, e.g. arg=Rocca1,Espirra
   -y      run year(s), e.g. arg=2003,2004-2005-2006
EOF

}


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

sitenames=()

for filename in ${met}*.nc; do

    sitename=$(echo $(basename -- "${filename}") |awk -F${suffix} '{print $1}')
    sitenames+=("${sitename}")

done

}


#######################################
# Default to all years
# Globals:
#   sitenames, runyrs
# Arguments:
#   None
# Returns:
#   array containing site run years
#######################################

all_years(){

years2run=()

for sitename in "${sitenames[@]}"; do # retrieve years

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
#   met, suffix, sitenames, runyrs
# Arguments:
#   None
# Returns:
#   sitenames, years2run,
#   (models, project)
#######################################

get_undeclared_variables(){

# NCPUs for parallelisation
if [[ -z ${NCPUs+x} ]]; then
    # unique physical cores
    unicore=$(grep ^cpu\\scores /proc/cpuinfo |uniq |awk '{print $4}')

    # check whether threads are physical or not 
    if [[ (( $(grep -E "cpu cores|siblings|physical id" /proc/cpuinfo \
               |xargs -n 11 echo |sort |uniq |awk '/^physical id/{print $7}' \
               |tail -1) -le $(grep -E "cpu cores|siblings|physical id" \
               /proc/cpuinfo |xargs -n 11 echo |sort |uniq \
               |awk '/^physical id/{print $11}' |tail -1) )) && \
          (( ${unicore} -lt $(nproc --all) )) ]]; then
        NCPUs=${unicore}

    else
        NCPUs=$((${unicore} / 2)) # hyperthreading, only take half the CPUs
    fi
fi

if [[ -z ${models+x} ]]; then # default is Control, ProfitMax
    models="default" 
fi

# is the project specified?
if [[ -z ${project+x} ]]; then
    project="None"
    runyrs="${data_dir}${sl}input${sl}fluxsites${sl}run_yrs_sites.txt"

else
    runyrs="${data_dir}${sl}input${sl}projects${sl}${project}${sl}"
    runyrs+="run_yrs_sites.txt"

    if [[ ! -f ${runyrs} ]]; then
        runyrs="${data_dir}${sl}input${sl}fluxsites${sl}run_yrs_sites.txt"
    fi
fi

# default to all sites in input/fluxsites/met
if [[ -z ${sitenames+x} ]]; then 
    all_sites
fi

# default to all years in input/fluxsites/run_yrs_sites.txt
if [[ ( ( -n ${years2run+x} ) && ( ${#years2run[@]} != ${#sitenames[@]} && \
          ${#years2run[@]} != 1 ) ) || ( -z ${years2run+x} ) ]]; then
    all_years
fi

# is the model type specified?
if [[ -z ${models+x} ]]; then
    models="default" 
fi

# is the photosynthesis model specified?
if [[ -z ${photomod+x} ]]; then
    photomod="default"
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

}


#######################################
# Get all combinations of site-year
# Globals:
#   sitenames, years2run
# Arguments:
#   None
# Returns:
#   combinations of site-year
#######################################

get_combinations(){

i=0 # initialise yearsofsite tracker
combinations=()

for sitename in "${sitenames[@]}"; do # loop over sites

    for yearsofsite in "${years2run[i]}"; do # loop over years of all sites

        # make string of years into array
        yearsofsite=$(echo ${yearsofsite} |sed 's,\-, ,g')
        yearsofsite=(${yearsofsite// / })

        for year2run in "${yearsofsite[@]}"; do # loop over years of this site

            combinations+=("${sitename}:${year2run}")

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
### The hydraulic stream's resolution is: ${resolution}
### The solving frequency is: ${fopt}
### The solving window is: ${window}
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

runfile="${this_dir}${sl}tmp${sl}irun${i}.txt"

if [[ ${1} == "spinup" ]]; then # make irunN.txt for spinup & run model
    ${mkrunfile} -s "${sitename}" -y "${year2run}" -m Control -P "${photomod}" \
                 -R "${project}" -o "${runfile}" >> ${log_file} 2>&1

elif [[ ${1} == "actual" ]]; then # make irunN.txt for actual & run model
    ${mkrunfile} -s "${sitename}" -a -y "${year2run}" -m "${models}" \
                 -r "${resolution}" -f "${fopt}" -w "${window}" \
                 -P "${photomod}" -R "${project}" -o "${runfile}" \
                 >> ${log_file} 2>&1
fi

if [[ ${1} == "spinup" ]]; then
    echo "Running spinup" ${sitename} ${year2run}
    echo "The ${sitename} ${year2run} ${photomod} spinup run" >> ${log_file}

elif [[ ${1} == "actual" ]]; then
    echo "Running actual" ${sitename} ${year2run}
    message="The ${sitename} ${year2run} ${photomod} ${resolution} ${fopt} "
    message+="${window} actual run"
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

to_rm=$(find ${this_dir}${sl}tmp${sl}"irun"*[0-9]*".txt" \
        |grep -v "${this_dir}${sl}tmp${sl}irun.txt" |grep -v ${log_file} |xargs)
rm ${to_rm}

}


#### execute functions #########################################################

# filter options that exist from those that don't
while getopts "hc:m:P:R:s:y:r:f:w:" OPTION; do

    case ${OPTION} in

        h)
            usage
            exit 1
            ;;

        c)
            NCPUs=${OPTARG}
            ;;

        m)
            models="${OPTARG}"
            ;;

        P)
            photomod=${OPTARG}
            ;;

        R)
            project=${OPTARG}
            ;;

        s)
            sitenames=()

            for i in $(echo ${OPTARG} |sed "s/,/ /g"); do

                sitenames+=("${i}")

            done

            ;;

        y)
            years2run=()

            for i in $(echo ${OPTARG} |sed "s/,/ /g"); do

                years2run+=("${i}")

            done

            if [[ ${#years2run[@]} != ${#sitenames[@]} ]]; then

                while [[ ${#years2run[@]} != ${#sitenames[@]} ]]; do

                    years2run+=("${years2run[-1]}")

                done

            fi

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

        ?)
            usage
            exit
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

# log file
log_file="${this_dir}${sl}tmp${sl}log.o$(date +%s)"

if [[ ! -d ${this_dir}${sl}tmp${sl} ]]; then
    mkdir ${this_dir}${sl}tmp
fi

# execute main & output run time in log
{ time main ; } >> ${log_file} 2>&1

