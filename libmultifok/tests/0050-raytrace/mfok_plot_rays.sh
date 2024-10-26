#! /bin/bash
# Last edited on 2024-10-22 03:54:23 by stolfi

rayFile="$1"; shift

frameEls=( ${rayfile/\// } )
title="${frameEls[2]}"

gnuplot <<EOF
  set term x11 size 800,800
  set size ratio -1
  set title "${title}"
  # Pixel center and half-side in scene coordinates (k is 0 for X, 1 for Y):
  ctrpix(k) = column(1+k)+radpix(k)
  radpix(k) = 0.5*column(3)
  # Ray hit point in scene coordinates (k is 0,1,2 for X,Y,Z):
  hit(k) = column(14+k)
  plot \
    "${rayFile}" using (ctrpix(0)):(ctrpix(1)):(radpix(0)):(radpix(1)) title "pixel" with boxxyerror lw 3 lc rgb '#775533', \
    "" using (hit(0)):(hit(1)) title "rays" with points pt 7 ps 0.75 lc rgb '#0022ff'
  pause mouse
EOF

  
  
