#! /bin/bash
# Last edited on 2014-06-14 02:54:47 by stolfilocal

name="$1"; shift  # Name of input and output file.

export GDFONTPATH=.

tmp=/tmp/$$

bimgf="${tmp}.png" # Raw image.
simgf="${name}-e.png" # Downscaled PNG file.

gnuplot <<EOF
set term png size 1600,800 font "arial,16"
set output "${bimgf}"
set size 1,1
set title "errors per bin"
set logscale y
set yrange [0.0000001:+0.200]

plot "${name}.txt" using 1:6 title "rms error" with linespoints pt 7 lc rgb '#005588'

quit
EOF

convert ${bimgf} -resize '50%' ${simgf}
rm -f ${bimgf}
display ${simgf}
