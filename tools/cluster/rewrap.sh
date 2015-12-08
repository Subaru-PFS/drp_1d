#!/bin/bash

#########################
# User-defined variables
#########################
amazed_version="15nov27"
#methods="correlationsolve chisquaresolve linematching linematching2 linemodel blindsolve fullsolve decisionaltree7 decisionaltreea decisionaltreeb"
methods="chisquaresolve"
dataset="batch6"
file_regex=".*_F\.fits"
error_file_regex=".*_ErrF\.fits"
################################################
# User should not need to change anything else.
################################################

run=`pwd | sed 's/.*run//'`
number_of_files=`ls -1 /data/cesam/rborges/$dataset/ | grep -E "${file_regex}" | wc -l`
max_index=$(($number_of_files-1))

variables="run dataset file_regex error_file_regex amazed_version number_of_files max_index"
echo "variables=\"$variables\"" > this_run.sh
for variable in $variables; do
    echo "$variable=${!variable}" >> this_run.sh
done
echo "methods=\"$methods\""   >> this_run.sh

echo "" > method_per_array.txt
for method in $methods; do
    job_name=`qsub -N run_$run_$method -o /net/fs/users/rborges/run$run/ -e /net/fs/users/rborges/run$run/ -q batch -v amazed_version=$amazed_version,run=$run,dataset=$dataset,file_regex=$file_regex,error_file_regex=$error_file_regex,method=$method /net/fs/users/rborges/run$run/run.sh -t 5`
    echo -n "$job_name" | sed 's/\[.*//' | tee -a method_per_array.txt
    echo " $method" | tee -a method_per_array.txt
done

echo -n "Waiting completion... ("
last_job_number=`echo $job_name | sed 's/\[.*//'`
last_job_file_name="run_$method.o$last_job_number-$max_index"
echo "$last_job_file_name)"

while [ ! -e $last_job_file_name ]; do
    sleep 10
done

echo "Last job finished!"
