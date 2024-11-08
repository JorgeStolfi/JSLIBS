#! /bin/bash
# Last edited on 2024-11-06 02:02:58 by stolfi

file="$1"; shift;

gnuplot <<EOF
set term X11 size 1200,800
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
pause mouse

EOF
