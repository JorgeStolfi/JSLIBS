#! /bin/bash
# Last edited on 2014-06-14 02:47:30 by stolfilocal

name="$1"; shift  # Name of input and output file.

export GDFONTPATH=.

tmp=/tmp/$$

bimgf="${tmp}.png" # Raw image.
simgf="${name}-v.png" # Downscaled PNG file.

gnuplot <<EOF
set term png size 1600,800 font "arial,16"
set output "${bimgf}"
set size 1,1
set title "Decoded value per encoded value"
set yrange [-4.500:+4.500]

plot "${name}.txt" using 1:2 notitle with linespoints pt 7 lc rgb '#008833'

quit
EOF

convert ${bimgf} -resize '50%' ${simgf}
rm -f ${bimgf}
display ${simgf}
