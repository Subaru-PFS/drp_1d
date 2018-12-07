#!/bin/bash

source this_run.sh
dataset_folder="/data/cesam/rborges/$dataset"

for method in $methods; do
    if [ "$methods_list" == "" ]; then
	methods_list="$method"
	durations_list="${method}_duration"
    else
	methods_list="$methods_list,$method"
	durations_list="$durations_list,${method}_duration"
    fi
done

rm {results,redshift}*.csv
echo "\"$run,$dataset,$templates_dir,$file_list,$error_file_list,$number_of_files,$amazed_version\",$number_of_files" > results.csv
echo "index,expected,$methods_list,$durations_list" >> results.csv

mv /data/cesam/rborges/$run/run*.{e,o}* .

for index in `seq 0 $max_index`; do
    prefix=`cat /data/cesam/rborges/$run/input_$index.spectrumlist | awk '{ print $1; }' | sed 's/\.fits//' | sed 's/_.*//'`
    expected=`grep -E "^$prefix[[:space:]]" $dataset_folder/Z-Mag.list | awk '{ print $2; }'`
    estimated_list=""
    duration_list=""
    for method in $methods; do
	duration_line=`grep -E "^real " run${run}_$method.e*-$index 2>&1`
	if [ $? -ne 0 ]; then
	    duration="err"
	else
	    duration=`echo $duration_line | tail -n1 | sed 's/.*real //'`
	fi
	if [ "$duration_list" == "" ]; then
	    duration_list="$duration"
	else
	    duration_list="$duration_list,$duration"
	fi
	estimated_line=`tail -n1 /data/cesam/rborges/$run/output_${index}_$method/redshift.csv 2>&1`
	if [ $? -eq 0 -a "$estimated_line" != "" ]; then
	    estimated=`echo $estimated_line | awk '{ print $2; }'`
	    tail -n1 /data/cesam/rborges/$run/output_${index}_$method/redshift.csv >> redshift_${method}.csv
	else
	    estimated="err"
	fi
	if [ "$estimated_list" == "" ]; then
	    estimated_list="$estimated"
	else
	    estimated_list="$estimated_list,$estimated"
	fi
    done
    echo "$index,$expected,$estimated_list,$duration_list" >> results.csv
done

mkdir torque_out
mv run${run}_*.{e,o}* torque_out/
mkdir amazed_out
mv /data/cesam/rborges/$run/* amazed_out/
rmdir /data/cesam/rborges/$run
