#! /bin/bash
# Last edited on 2023-03-29 19:36:23 by stolfi

usage="$0 PROBLEM INTEGRATOR SECONDS"

prbl="$1"; shift
intg="$1"; shift
secs="$1"; shift;

gnuplot <<EOF
set terminal X11
plot \
  "out/${prbl}-${intg}-merr.plot" using 1:2 with points, \
  "out/${prbl}-${intg}-eerr.plot" using 1:2 with points
pause $secs
quit
EOF

exit 0
