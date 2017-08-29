#!/bin/bash

gnuplot <<- EOF
  set xlabel "Number of threads"
  set ytics  nomirror
  set ylabel "Time (seconds)"
  set yrange[1:]
  set y2range[0:]
  set xrange[1:32]
  set grid
  set y2tics  nomirror
  set y2label "Speedup (factor)"
  set key top center
  set terminal pngcairo dashed size 800,400 enhanced
  set output 'scalability.png'
  set style line 1 lt 1 lw 1 pt 1 linecolor rgb "red"
  set style line 2 lt 1 lw 1 pt 1 linecolor rgb "navy"
  set style line 3 lt 2 lw 1 pt 1 linecolor rgb "navy"
  plot "scalability.dat" using 1:2 title 'Execution Time' with linespoints ls 1 axes x1y1, "scalability.dat" using 1:(31.2136/\$2) title 'Speedup' with linespoints ls 2 axes x2y2, "scalability.dat" using 1:1 title 'Perfect speedup' with linespoints ls 3 axes x2y2
EOF

