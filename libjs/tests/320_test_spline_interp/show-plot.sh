#! /bin/bash
# Last edited on 2019-04-09 19:33:05 by jstolfi

file="$1"; shift;
if [[ -s ${file} ]]; then 

  gnuplot <<EOF
  set term postscript eps color "TimesRoman" 20
  set size 3,1
  set title "${file}"
  # Plot the interpolation kernel:
  set output "${file%.*}-ker.eps"
  set xrange [4:16]
  plot \
    "${file}" using 1:2 notitle with linespoints lt 1 pt 7 ps 0.5, \
    "${file}" using 1:3 notitle with points pt 7 lt 0
  
  # Plot the interpolation kernel:
  set output "${file%.*}-tst.eps"
  set xrange [0:*]
  plot \
    "${file}" using 1:2 notitle with linespoints lt 1 pt 7 ps 0.1, \
    "${file}" using 1:3 notitle with points pt 7 lt 0 ps 0.4
EOF

atril "${file%.*}-ker.eps"
atril "${file%.*}-tst.eps"

fi
