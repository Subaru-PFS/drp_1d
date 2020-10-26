#!/bin/bash

get_info() {
    tmp=/tmp/monitor_torque.sh.tmp
    showq > $tmp
    jobs_in_flight=`cat $tmp | grep rborges | wc -l`
    num_jobs_held=`cat $tmp | grep rborges | grep Hold | wc -l`
    num_jobs_deferred=`cat $tmp | grep rborges | grep Deferred | wc -l`
    jobs_idle=`cat $tmp | grep rborges | grep Idle | wc -l`
    jobs_running=`cat $tmp | grep rborges | grep Running | wc -l`
}

print_info() {
    now=`date +"%F_%T"`
    printf "%s \t %d \t %d \t %d \t\t %d \t %d \n" $now $jobs_in_flight $num_jobs_held $num_jobs_deferred $jobs_idle $jobs_running
}

print_header() {
    now=`date +"%F_%T"`
    printf "$now \t Total \t Hold \t Deferred \t Idle \t Running\n"
    for i in {0..71}; do
	echo -n "="
    done
    echo
}

count=0
while [ true ]; do
    modulo=$(($count%35))
    if [ $modulo -eq 0 ]; then
	print_header
    fi
    get_info
    print_info
    sleep 60
    count=$(($count+1))
done
