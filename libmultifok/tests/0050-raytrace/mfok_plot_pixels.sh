#! /bin/bash
# Last edited on 2024-12-15 22:40:13 by stolfi

pixelFile="$1"; shift

pixel="${pixelFile##*/}"
title="${pixel}" 

gnuplot <<EOF
  set term x11 noenhanced size 1200,800
  set title "${title}"

  iPix(k) = column(1+k)
  pixSize(k) = column(3)
  pCtr(k) = column(4+k)
  zFoc(k) = column(7)
  zDep(k) = column(8)
  shrp(k) = column(9)
  vBlr(k) = 1/shrp(k)
  hAvg(k) = column(10)
  hDev(k) = column(11)
  plot \
    "${pixelFile}" using (zFoc(0)):(shrp(1)) title "shrp" with linespoints pt 7 ps 1.50 lc rgb '#ff4400', \
    ""             using (zFoc(0)):(hAvg(1)) title "hAvg" with linespoints pt 7 ps 0.75 lc rgb '#008800', \
    ""             using (zFoc(0)):(hDev(1)) title "hDev" with linespoints pt 7 ps 0.75 lc rgb '#0022ff'
  pause mouse
EOF

  
  
