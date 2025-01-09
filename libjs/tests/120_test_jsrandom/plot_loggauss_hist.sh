#! /bin/bash
# Last edited on 2025-01-02 06:32:59 by stolfi

# Usage: "plot-hist.sh {FILE}"
# Plots histogram {FILE}.his to {FILE}.png

name="$1"; shift

hfile="out/${name}.his"
pfile="out/${name}.png"

if [[ ! ( -s ${hfile} ) ]]; then echo "** ${hfile} not found" 1>&2; exit 1; fi
 
gnuplot <<EOF
set term png small truecolor background rgb '#ffffff'
set output "${pfile}"
set xrange [-0.25:+37.50]
set yrange [-0.10:+0.40]

# Average and deviation of log-gauss var from test program:
avg = 7.0
dev = 2.0
# Standard parameters of log-normal distr:
pi = 3.1415926
rr = dev/avg
mu = log(avg/sqrt(1.0 + rr*rr))
si = sqrt(log(1.0 + rr*rr))
# PDF of log-normal distr:
loggauss(x) = exp(-(log(x) - mu)**2/(2*si**2))/(x*si*sqrt(2*pi))
# Center of histogram bin:
avgc(i,j) = (column(i) + column(j))/2
plot \
  "${hfile}" using (avgc(3,4)):(loggauss(avgc(3,4))) title "expected" with lines lt 1 lc rgb '#ff0000', \
  "${hfile}" using 3:5 title "observed" with histeps lt 1 lc rgb '#008800'
EOF
display "${pfile}"
