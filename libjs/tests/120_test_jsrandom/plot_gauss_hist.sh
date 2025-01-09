#! /bin/bash
# Last edited on 2025-01-02 06:32:33 by stolfi

# Usage: "plot-hist.sh {FILE}"
# Plots histogram {FILE}.his to {FILE}.png

name="$1"; shift

hfile="out/${name}.his"
pfile="out/${name}.png"

if [[ ! ( -s ${hfile} ) ]]; then echo "** ${hfile} not found" 1>&2; exit 1; fi

gnuplot <<EOF
set term png small truecolor background rgb '#ffffff'
set output "${pfile}"
set xrange [-6.25:+6.25]
set yrange [-0.20:+0.60]
avgc(i,j) = (column(i) + column(j))/2
gauss(x) = exp(-x*x/2)/sqrt(2*3.1415926)
plot \
  "${hfile}" using (avgc(3,4)):(gauss(avgc(3,4))) title "expected" with lines lt 1 lc rgb '#ff0000', \
  "${hfile}" using 3:5 title "observed" with histeps lt 1 lc rgb '#008800'
EOF
display "${pfile}"
