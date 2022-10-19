#! /bin/bash
# Last edited on 2021-12-19 23:25:51 by stolfi

datafile="$1"; shift

plotfile="${datafile/.dat/.png}"

tplotfile="out/.tmp.png"

gnuplot <<EOF
  set term pngcairo size 1200,1200
  set output "${tplotfile}"
  set xrange [-5:+5]
  set yrange [-5:+5]
  set zrange [-0.001:]
  set hidden3d
  splot \
    "${datafile}" using 1:2:3 with points pt 7 ps 2, \
    "${datafile}" using 1:2:(0) with points pt 6 ps 2 
  
EOF

convert out/.tmp.png -resize '50%' out/f2.png 
display out/f2.png 
