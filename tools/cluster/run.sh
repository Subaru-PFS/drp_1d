#!/bin/bash

#PBS -l walltime=30:00

mkdir -p /data/cesam/rborges/$run/output_$PBS_ARRAYID

file_name=`head -n $(($PBS_ARRAYID+1)) $file_list | tail -n 1`

if [ "$error_file_list" != "" ]; then
    error_file_name=`head -n $(($PBS_ARRAYID+1)) $error_file_list | tail -n 1`
    echo -e "$file_name	$error_file_name" > /data/cesam/rborges/$run/input_$PBS_ARRAYID.spectrumlist
else
    echo -e "$file_name" > /data/cesam/rborges/$run/input_$PBS_ARRAYID.spectrumlist
fi

random_sleep() {
    sleep $((10+$RANDOM%10))
}

# scratch files
###############
scratch=/tmp/amazed
mkdir $scratch
job_input=$scratch/input_$PBS_ARRAYID.spectrumlist
random_sleep
cp /data/cesam/rborges/$run/input_$PBS_ARRAYID.spectrumlist $job_input
job_templates_dir=$scratch/$templates_dir
if [ ! -e $job_templates_dir ]; then
    random_sleep
    cp -r /data/cesam/rborges/$templates_dir $job_templates_dir
fi
job_spectrum=$scratch/$dataset
if [ ! -e $job_spectrum ]; then
    cp -r /data/cesam/rborges/$dataset $job_spectrum
fi
job_linecatalog=$scratch/linecatalogamazedair_B5.txt
if [ ! -e $job_linecatalog ]; then
    cp /data/cesam/rborges/linecatalogamazedair_B5.txt $job_linecatalog
fi
job_parameters_file=$scratch/$parameters_file
if [ ! -e $job_parameters_file ]; then
    cp /data/cesam/rborges/$parameters_file $job_parameters_file
fi

# amazed
time -p /data/cesam/rborges/amazed_$amazed_version --input $job_input --templatedir $job_templates_dir --spectrum $job_spectrum -y $job_linecatalog -o /data/cesam/rborges/$run/output_${PBS_ARRAYID}_$method -m $method --parameters $job_parameters_file --thread-count 1
