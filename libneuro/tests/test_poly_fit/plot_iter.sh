#! /bin/bash
# Last edited on 2013-11-23 00:33:08 by stolfilocal

initf="$1"; shift
iterf="$1"; shift
echo "== ${iterf}"

export GDFONTPATH=.

tmp=/tmp/$$

bimgf="${tmp}.png" # Raw image.
simgf="${iterf%%.*}.png" # Downscaled PNG file.

gnuplot <<EOF
set term png size 1600,1600 font "arial,16"
set output "${bimgf}"
set size 1,1
set title "${iterf}"
set yrange [-2.10:+2.10]

plot "${initf}" using 2:3 title "true" with lines lc rgb '#0000ff', \
     "${initf}" using 2:4 title "data" with points ps 2.0 pt 6 lc rgb '#005588', \
     "${iterf}" using 2:4 title "cook" with points ps 1.5 pt 7 lc rgb '#ff0000', \
     "${iterf}" using 2:3 title "appr" with lines lc rgb '#ff00ff'

quit
EOF

convert ${bimgf} -resize '50%' ${simgf}
rm -f ${bimgf}
