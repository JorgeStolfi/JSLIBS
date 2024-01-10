#! /bin/bash
# Last edited on 2013-12-17 01:03:39 by stolfilocal

# Usage: "plot-hist.sh {FILE}"
# Plots {FILE}.his to {FILE}.png

name="$1"; shift

hfile="out/${name}.his"
pfile="out/${name}.png"

gnuplot -background white <<EOF
set term png small truecolor
set output "${pfile}"
set xrange [-0.75:]
set yrange [-0.20:+1.20]
const(v,z) = v
plot \
  const(1,x) with lines lt 1 lc rgb '#ff0000', \
  "${hfile}" using 1:3 with histeps lt 1 lc rgb '#008800'
EOF
display "${pfile}"
