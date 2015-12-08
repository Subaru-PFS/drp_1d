#!/bin/bash

range_min=$1
range_max=$2
suffix='_F'

if [ "$range_min" == "" -o \
    "$range_max" == "" ]; then
    echo "Filter input.spectrumlist into a filtered.spectrumlist that is within the input redshift interval. Interval is [range_min, range_max]. Requires that input.spectrumlist and Z-Mag.list be on the current working directory. Supposes input.spectrumlist lines are of the format <anything>$suffix.fits."
    echo "USAGE:"
    echo "$0 <range_min> <range_max>"
    exit 1
fi

if [ -e filtered.spectrumlist ]; then
    rm filtered.spectrumlist
fi

while read line; do
    file_name=`echo "$line" | awk '{ print $1; }'`
    index=`echo $file_name | sed "s/$suffix\.fits.*//"` # <index>_F.fits
    #echo $index
    expected_redshift=`grep -E "^$index " Z-Mag.list | awk '{ print $2; }'` 
    #echo $expected_redshift
    if [ "$expected_redshift" == "" ]; then
	echo "Z-Mag mising $index"
        continue
    fi
    min_test=`echo "$expected_redshift >= $range_min" | bc`
    if [ $min_test -eq 0 ]; then
	continue
    fi
    max_test=`echo "$expected_redshift <= $range_max" | bc`
    if [ $max_test -eq 0 ]; then
	continue
    fi
    # At this point, redshift is within interval.
    echo "$line" >> filtered.spectrumlist    
done < input.spectrumlist
