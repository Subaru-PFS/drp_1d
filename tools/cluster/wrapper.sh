#!/bin/bash

#########################
# User-defined variables
#########################
amazed_version="151126"
complex_methods="blindsolve fullsolve decisionaltree7 decisionaltreea decisionaltreeb"
operators="correlationsolve chisquare2solve linematching linematching2 linemodel"
dataset="batch6"
file_list="/data/cesam/rborges/batch6/batch6_F_file_list.txt"
error_file_list="/data/cesam/rborges/batch6/batch6_ErrF_file_list.txt"
templates_dir="ExtendedGalaxyEL3"
parameters_file="all_methods_03.json"
################################################
# User should not need to change anything else.
################################################

# Sanity
########
if [ ! -e $file_list ]; then
    echo "$file_list does not exist!"
    exit 1
fi
if [ ! -e $error_file_list ]; then
    echo "$error_file_list does not exist!"
    exit 1
fi
if [ ! -e /appli/cesam/rborges/amazed_$amazed_version ]; then
    echo "/appli/cesam/rborges/amazed_$amazed_version does not exist!"
    exit 1
fi
if [ ! -e /data/cesam/rborges/$parameters_file ]; then
    echo "/data/cesam/rborges/$parameters_file does not exist!"
    exit 1
fi

methods="$operators $complex_methods"
burst=31
run=`pwd | sed 's/.*run//'`
number_of_files=`cat ${file_list} | wc -l`
max_index=$(($number_of_files-1))

variables="run dataset file_list error_file_list amazed_version number_of_files max_index templates_dir parameters_file"
echo "variables=\"$variables\"" > this_run.sh
for variable in $variables; do
    echo "$variable=${!variable}" >> this_run.sh
done
echo "methods=\"$methods\""   >> this_run.sh

mkdir -p /data/cesam/rborges/${run}

echo "" > method_per_array.txt
for method in $methods; do

    final_burst=1
    burst_max_index=-1
    while [ $final_burst -ne 0 ]; do
	burst_min_index=$(($burst_max_index+1))
	burst_max_index=$(($burst_min_index+$burst))
	if [ $burst_max_index -ge $max_index ]; then
	    final_burst=0
	    burst_max_index=$max_index
	fi
	job_name=`qsub -N run${run}_${method} -o /data/cesam/rborges/${run} -e /data/cesam/rborges/${run} -q batch -v amazed_version=$amazed_version,run=$run,dataset=$dataset,file_list=$file_list,error_file_list=$error_file_list,method=$method,templates_dir=$templates_dir,parameters_file=$parameters_file /net/fs/users/rborges/run$run/run.sh -t $burst_min_index-$burst_max_index`
	echo -n "$job_name" | sed 's/\[.*//' | tee -a method_per_array.txt
	echo " $method $burst_min_index-$burst_max_index" | tee -a method_per_array.txt

	echo -n "Waiting completion... ("
	burst_half_index=$((($burst_max_index+$burst_min_index)/2))
	half_job_number=`echo $job_name | sed 's/\[.*//'`
	half_job_file_name="run${run}_$method.o$half_job_number-$burst_half_index"
	echo "$half_job_file_name)"
	while [ ! -e /data/cesam/rborges/$run/$half_job_file_name ]; do
	    sleep 10
	done
    done # burst

done # method
echo -n "Waiting completion... ("
last_job_file_name="run${run}_$method.o$half_job_number-$max_index"
echo "$last_job_file_name)"
while [ ! -e /data/cesam/rborges/$run/$last_job_file_name ]; do
    sleep 10
done
echo "Everything finished!"
