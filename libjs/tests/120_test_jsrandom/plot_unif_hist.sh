#! /bin/bash
# Last edited on 2020-12-06 17:51:05 by jstolfi

# Usage: "plot-hist.sh {FILE}"
# Plots histogram {FILE}.his to {FILE}.png

name="$1"; shift

hfile="out/${name}.his"
pfile="out/${name}.png"

gnuplot -background white <<EOF
set term png small truecolor
set output "${pfile}"
set xrange [-0.20:+1.20]
set yrange [-0.20:+1.20]
set xzeroaxis
const(v,z) = ((z < -0.01) || (z > 1.001) ? 0/0 : v)
plot \
  const(1,x) title "expected" with lines lt 1 lc rgb '#ff0000', \
  "${hfile}" using 3:5 title "observed" with histeps lt 1 lc rgb '#008800'
EOF
display "${pfile}"
