#!/usr/bin/env bash

if [[ $# != 3 ]] ; then
    echo "Script usage:"
    echo "./reductions.sh binary_name input_file num_threads"
    exit -1
fi

BIN=$1
input=$2
NUM_TH=$3

tmpf="tmptest.pgm"
output="output.txt"

NUM_REP=5
res=0

for it in 1 2 3 4 5 10 15 20 25 50 100 150 200 20000
do
    #echo "Testing $inputs with $it iterations between reductions:"
    echo " " > $output
    for re in $(seq $NUM_REP)
    do
        mpirun -n $NUM_TH ./$BIN $input 10000 $it $tmpf  >> $output
    done
    # Get the average of iteration time
    res=$(cat $output | grep "Executed" | awk 'BEGIN{acc=0;n=0}{acc+=$5;n+=1}END{print acc/n}')
    echo $it $res
done

rm $tmpf output.txt
