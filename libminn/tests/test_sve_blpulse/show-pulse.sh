#! /bin/bash
# Last edited on 2009-02-22 19:03:00 by stolfi

# Plots the bandlimited pulse by {test_sve_bpulse}

file="$1"; shift;
if [[ -r ${file} ]]; then
  gnuplot << EOF
  set terminal X11
  plot \
    "${file}" using 1:2 with linespoints
  pause 300
EOF
else
  echo "file ${file} missing" 1>&2
fi

