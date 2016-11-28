#!/usr/bin/env bash

if [[ $# != 3 ]] ; then
    echo "Script usage:"
    echo "./image_size.sh binary_name num_threads inputs_root"
    exit -1
fi

BIN=$1
NUM_TH=$2
ir=$3

tmpf="tmptest.pgm"
output="output.txt"

NUM_REP=5
res=0
red=0
wri=0

for s in 100 300 500 700 900 1100 1300 1500 1700 1900 2100 2300 2500
do
    echo " " > $output
    for re in $(seq $NUM_REP)
    do
        mpiexec -n $NUM_TH ./$BIN $ir${s}x${s}.pgm 200 200 $tmpf  >> $output
    done
    # Get the average of iteration time
    red=$(cat $output | grep "Data readed" | awk 'BEGIN{acc=0;n=0}{acc+=$6;n+=1}END{print acc/n}')
    wri=$(cat $output | grep "Data gathered" | awk 'BEGIN{acc=0;n=0}{acc+=$6;n+=1}END{print acc/n}')
    res=$(cat $output | grep "Executed" | awk 'BEGIN{acc=0;n=0}{acc+=$5;n+=1}END{print acc/n}')
    echo $s $res $red $wri
done

rm $tmpf $output
