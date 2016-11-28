#!/usr/bin/env bash

if [[ $# != 3 ]] ; then
    echo "Script usage:"
    echo "./ctest.sh binary_name path_inputs path_correct_results"
    exit -1
fi


BIN=$1
PI=$2
PR=$3
tmpf="tmptest.pgm"
out="out.txt"
ave="averages.txt"
inputs=('edge192x128.pgm' 'edge768x768.pgm' )
for i in ${inputs[@]}
do
   for t in 1 2 3 4 8
   do
     echo "Testing $i with $t threads:"
      mpirun -n $t ./$BIN $PI/$i 0 100 $tmpf > $out
      cat $out | grep 'average' > $ave
      cmp $tmpf $PR/$i
      cmp $ave $PR/$i.averages
   done
done

inputs=('edge192x128.pgm' )
for i in ${inputs[@]}
do
   t=9
   echo "Testing $i with $t threads:"
   eval mpirun -n $t ./$BIN $PI/$i 0 5 $tmpf &> /dev/null
   ret_code=$?
   if [ $ret_code != 0 ]; then
     echo "Test ok!"   
   else
     echo "Test failed! Image should not be divisible in 9 processes"
   fi
done

rm $tmpf $out $ave

