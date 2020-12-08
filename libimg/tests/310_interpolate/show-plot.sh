#! /bin/bash
# Last edited on 2012-01-27 04:02:45 by stolfilocal

file="$1"; shift;

gnuplot -geometry '1200x800' <<EOF
set term X11
set size 1,1
set surface
set view 80,30
set hidden3d
set xyplane at -2
set title "${file}"
splot \
  "${file}" using 1:2:3                 title "order = -1" with lines lt 2, \
  "${file}" using 1:2:(column(4)+2)     title "order =  0" with lines lt 1, \
  "${file}" using 1:2:(column(5)+4)     title "order =  1" with lines lt 3, \
  "${file}" using 1:2:(0.1*column(6)-2) title "domain" with lines lt 4
pause 30

EOF
