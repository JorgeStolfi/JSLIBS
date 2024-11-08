#! /bin/bash
# Last edited on 2024-11-05 11:21:57 by stolfi

# Usage: "plot-hist.sh {FILE}"
# Plots histogram {FILE}.his to {FILE}.png

name="$1"; shift

hfile="out/${name}.his"
pfile="out/${name}.png"

gnuplot<<EOF
set term png small truecolor background rgb '#ffffff'
set output "${pfile}"
set xrange [-0.75:]
set yrange [-0.20:+1.20]
const(v,z) = v
plot \
  const(1,x) with lines lt 1 lc rgb '#ff0000', \
  "${hfile}" using 1:3 with histeps lt 1 lc rgb '#008800'
EOF
display "${pfile}"
