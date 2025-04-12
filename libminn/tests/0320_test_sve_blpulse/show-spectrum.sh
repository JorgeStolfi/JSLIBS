#! /bin/bash
# Last edited on 2025-03-19 15:11:12 by stolfi

# Plots the power spectrum found by {test_sve_bpulse}

file="$1"; shift;

if [[ -r ${file} ]]; then
  gnuplot << EOF
  set terminal X11
  set logscale y
  set title "spectrum (${file})"
  plot \
    "${file}" using 1:2 title "pwr" with linespoints, \
    "${file}" using 1:3 title "bad" with linespoints
  pause mouse
EOF
else
  echo "file ${file} missing" 1>&2
fi
