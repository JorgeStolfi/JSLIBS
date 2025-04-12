#! /bin/bash
# Last edited on 2025-04-02 18:36:39 by stolfi

dfile="$1"; shift

export GDFONTPATH="${HOME}/ttf"

gnuplot <<EOF

set term X11

set xlabel "distance (px)"
set ylabel "energy"

set yrange [1.0e-10:]
set xzeroaxis ls 1 lc rgb '#ffaa55'
set logscale y

plot \
  "${dfile}" using 1:2 title "energy"    with linespoints pt 7 ps 0.75 lc rgb '#ff0000', \
  "${dfile}" using 1:3 title "rterm"     with linespoints pt 7 ps 0.75 lc rgb '#009900', \
  "${dfile}" using 1:4 title "eterm"     with linespoints pt 7 ps 0.75 lc rgb '#aa00ff'

pause mouse

EOF
