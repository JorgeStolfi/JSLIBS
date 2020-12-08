#! /bin/bash
# Last edited on 2014-06-14 02:55:43 by stolfilocal

name="$1"; shift  # Name of input and output file.

export GDFONTPATH=.

tmp=/tmp/$$

bimgf="${tmp}.png" # Raw image.
simgf="${name}-p.png" # Downscaled PNG file.

gnuplot <<EOF
set term png size 1600,800 font "arial,16"
set output "${bimgf}"
set size 1,1
set title "Distribution of Gaussian samples per encoded value"
set yrange [-0.0000001:]

plot "${name}.txt" using 1:4 notitle with linespoints pt 7 lc rgb '#008833'

quit
EOF

convert ${bimgf} -resize '50%' ${simgf}
rm -f ${bimgf}
display ${simgf}
