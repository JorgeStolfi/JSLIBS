#! /bin/bash
# Last edited on 2013-10-26 01:48:18 by stolfilocal

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

gv "${file%.*}-ker.eps"
gv "${file%.*}-tst.eps"

fi
