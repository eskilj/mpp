#!/bin/bash

gnuplot <<- EOF
  set xlabel "Iterations between reductions"
  set ylabel "Time (seconds)"
  set yrange[1:]
  set xrange[1:]
  set grid
  set key top right
  set terminal pngcairo dashed size 800,400 enhanced
  set output 'reductions.png'
  set style line 1 lt 1 lw 1 pt 1 linecolor rgb "navy"
  set style line 2 lt 2 lw 1 pt 0 linecolor rgb "red"
  plot "reductions.dat" using 1:2 title 'Execution Time' with linespoints ls 1, "reductions.dat" using 1:3 title 'Execution time without reductions' with linespoints ls 2
EOF

