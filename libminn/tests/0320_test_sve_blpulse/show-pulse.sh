#! /bin/bash
# Last edited on 2025-03-19 15:11:02 by stolfi

# Plots the bandlimited pulse by {test_sve_bpulse}

file="$1"; shift;
if [[ -r ${file} ]]; then
  gnuplot << EOF
  set terminal X11
  set title "amplitude (${file})"
  plot \
    "${file}" using 1:2 title "amp" with linespoints
  pause mouse
EOF
else
  echo "file ${file} missing" 1>&2
fi

