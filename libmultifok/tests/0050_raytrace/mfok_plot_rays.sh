#! /bin/bash
# Last edited on 2024-12-15 22:32:56 by stolfi

rayFile="$1"; shift

pixel="${rayFile##*/}"
frame="${rayFile%/*}"; frame="${frame##*/}"
title="${frame}/${pixel}" 

gnuplot <<EOF
  set term x11 noenhanced size 800,800
  set size ratio -1
  set title "${title}"
  # Sample point {pSmp} in image coordinates: (k is 0 for X, 1 for Y, 2 for Z):
  iPix(k) = column(1+k)
  iSmp(k) = column(4+k)
  step(k) = column(6)
  wSmp(k) = column(7)
  # Sample point {pSmp} in scene coordinates: (k is 0 for X, 1 for Y, 2 for Z):
  sizePix(k) = column(3)
  radPix(k) = 0.5*sizePix(k)
  ctrPix(k) = (iPix(k)+0.5)*sizePix(k)
  pSmp(k) = ctrPix(k) + step(k)*iSmp(k)*sizePix(k)
  # Pixel center and half-side in scene coordinates (k is 0 for X, 1 for Y, 2 for Z):
  pRay(k) = column(8+k)
  dRay(k) = column(11+k)
  wRay(k) = column(14)
  pHit(k) = column(15+k)
  hHit(k) = column(18)
  vBlr(k) = column(19)
  plot \
    "${rayFile}" using (ctrPix(0)):(ctrPix(1)):(radPix(0)):(radPix(1)) title "pixel" with boxxyerror lw 3 lc rgb '#775533', \
    "" using (pRay(0)):(pRay(1)) title "sampoints" with points pt 7 ps 1.50 lc rgb '#ff4400', \
    "" using (pHit(0)):(pHit(1)) title "rays" with points pt 7 ps 0.75 lc rgb '#0022ff'
  pause mouse
EOF

  
  
