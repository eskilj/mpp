#!/bin/bash

gnuplot <<- EOF
set xlabel "Size of the squared images (pixels)"
  set ylabel "Time (seconds)"
  set key top right
  set grid
  set terminal pngcairo dashed size 800,400 enhanced
  set output 'imagesize.png'
  set style line 1 lt 1 lw 1 pt 1 linecolor rgb "navy"
  set style line 2 lt 1 lw 1 pt 1 linecolor rgb "red"
  set style line 3 lt 1 lw 1 pt 1 linecolor rgb "green"
  plot "imagesize.dat" using 1:2 title 'Read data and distribute' with linespoints ls 1, "imagesize.dat" using 1:3 title 'Gather data and write' with linespoints ls 2, "imagesize.dat" using 1:4 title 'Reverse edge-detection algorithm' with linespoints ls 3

EOF

