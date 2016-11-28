#!/bin/bash

gnuplot <<- EOF
  set xlabel "Number of threads"
  set ylabel "Parallel efficiency"
  set xrange[1:32]
  set yrange[0:2]
  set grid
  set key top center
  set terminal pngcairo dashed size 800,400 enhanced
  set output 'pefficiency.png'
  set style line 1 lt 1 lw 1 pt 1 linecolor rgb "red"
  plot "scalability.dat" using 1:(31.2136/(\$1*\$2)) title 'Parallel efficiency' with linespoints ls 1 
EOF

