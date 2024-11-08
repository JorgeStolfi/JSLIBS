#! /bin/bash
# Last edited on 2024-11-05 11:21:28 by stolfi

# Usage: "plot-hist.sh {FILE}"
# Plots histogram {FILE}.his to {FILE}.png

name="$1"; shift

hfile="out/${name}.his"
pfile="out/${name}.png"

gnuplot <<EOF
set term png small truecolor background rgb '#ffffff'
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
