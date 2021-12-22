#! /bin/bash
# Last edited on 2012-01-21 01:49:38 by stolfi

gnuplot <<EOF
  set term png medium size 1200,1200
  set output "out/.tmp.png"
  set hidden3d
  splot "out/f2.dat" using 1:2:3 with lines 
EOF

convert out/.tmp.png -resize '50%' out/f2.png 
display out/f2.png 
