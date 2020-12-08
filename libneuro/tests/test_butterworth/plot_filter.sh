#! /bin/bash 
# Last edited on 2013-06-14 08:59:39 by stolfilocal

datfile="$1"; shift;

export GDFONTPATH=.

gnuplot <<EOF

set term x11
set size 1,0.80
Gmin = 1.0e-5
Gmax = 1
set logscale y
fmin = 0.01
fmax = 100
set logscale x
set grid xtics
set grid ytics
set xrange [(0.99*fmin):(1.01*fmax)]
set yrange [(0.99*Gmin):(1.01*Gmax)]
plot \
  "${datfile}" using 1:(column( 2)**2) title 'n=1' with lines, \
  ""           using 1:(column( 3)**2) title 'n=2' with lines, \
  ""           using 1:(column( 4)**2) title 'n=3' with lines, \
  ""           using 1:(column( 5)**2) title 'n=4' with lines, \
  ""           using 1:(column( 6)**2) title 'n=5' with lines, \
  ""           using 1:(column( 7)**2) title 'n=6' with lines, \
  ""           using 1:(column( 8)**2) title 'n=7' with lines, \
  ""           using 1:(column( 9)**2) title 'n=8' with lines, \
  ""           using 1:(column(10)**2) title 'n=9' with lines
pause 300
EOF
