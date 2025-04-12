#! /bin/bash
# Last edited on 2025-03-19 15:15:58 by stolfi

fname="$1"; shift

gnuplot <<EOF
set term x11
plot \
  "${fname}" using 2:3  with linespoints pt 7 lt 1 lc rgb '#ff0000', \
  ""         using 2:4  with linespoints pt 7 lt 1 lc rgb '#aa5500', \
  ""         using 2:5  with linespoints pt 7 lt 1 lc rgb '#226600', \
  ""         using 2:6  with linespoints pt 7 lt 1 lc rgb '#00aa00', \
  ""         using 2:7  with linespoints pt 7 lt 1 lc rgb '#008888', \
  ""         using 2:8  with linespoints pt 7 lt 1 lc rgb '#0033ff', \
  ""         using 2:9  with linespoints pt 7 lt 1 lc rgb '#7700ff', \
  ""         using 2:10 with linespoints pt 7 lt 1 lc rgb '#aa0044'
pause mouse
quit
EOF
