#!/usr/bin/env bash

if [[ $# != 2 ]] ; then
    echo "Script usage:"
    echo "./scalability.sh binary_name input_file"
    exit -1
fi

BIN=$1
input=$2

tmpf="tmptest.pgm"
output="output.txt"

NUM_REP=5
res=0

for t in 1 2 4 8 16 32 64
do
    echo " " > $output
    for re in $(seq $NUM_REP)
    do
        mpiexec -n $t ./$BIN $input 20000 20 $tmpf >> $output
    done
    # Get the average of iteration time
    res=$(cat $output | grep "Executed" | awk 'BEGIN{acc=0;n=0}{acc+=$5;n+=1}END{print acc/n}')
    echo $t $res
done

rm $tmpf output.txt
