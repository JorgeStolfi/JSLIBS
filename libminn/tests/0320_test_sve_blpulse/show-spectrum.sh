#! /bin/bash
# Last edited on 2012-01-23 04:30:18 by stolfilocal

# Plots the power spectrum found by {test_sve_bpulse}

file="$1"; shift;

if [[ -r ${file} ]]; then
  gnuplot << EOF
  set terminal X11
  set logscale y
  plot \
    "${file}" using 1:2 with linespoints, \
    "${file}" using 1:3 with linespoints
  pause 300
EOF
else
  echo "file ${file} missing" 1>&2
fi
