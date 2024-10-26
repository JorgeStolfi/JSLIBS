#! /bin/bash 
# Last edited on 2023-01-22 13:26:26 by stolfi

rawFile="$1"; shift # File with ray tilt data.

tmp="/tmp/$$"

tempFile="${tmp}-d.txt"

cat ${rawFile} \
  | sed \
      -e 's:tilt\[[ 0-9]*\] *= *[()]::g' \
      -e 's:[()] *wr *=::g' \
  > ${tempFile}

export GDFONTPATH="${HOME}/tt-fonts"

gnuplot << EOF
set term X11 size 600,600 
set size ratio -1 square

wr(k) = 0.1*column(k)

plot "${tempFile}" using 1:2         title "positions" with points pt 7 ps 0.75 lc rgb '#008800', \
     ""            using 1:2:(wr(3)) title "weights"   with circles             lc rgb '#882200'

pause mouse
quit
EOF
