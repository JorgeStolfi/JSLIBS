#! /bin/bash
# Last edited on 2024-12-20 13:28:48 by stolfi

usage="${0} '{KIND}' '{TITLE1}' '{TITLE2}' {TABLE_NAME}"

# Where {KIND} is {sRGB} or {BT.709}.

kind="$1"; shift
escale="$1"; shift
dscale="$1"; shift
title1="$1"; shift
title2="$1"; shift
tname="$1"; shift

infile="${tname}.txt"
epsfile="${tname}.eps"

gnuplot <<EOF
set term postscript eps color solid linewidth 4.0 "Courier-Bold" 30
set output "${epsfile}"
set size 2.85,3
set size ratio -1
set xrange [:+2.00]
set logscale x
set yrange [-0.05:+1.05]
set title "${title1}\n${title2}"
set xzeroaxis
set yzeroaxis
set key outside
escale = ${escale}
dscale = ${dscale}
# test terminal
plot \
  "${infile}" using 1:2 title "${kind}" with lines lc rgb '#ff0000', \
  "${infile}" using 1:3 title "generic" with lines lc rgb '#0055ff', \
  "${infile}" using 1:(0.5 + escale*column(4)) title "diff x ${escale}" with lines lc rgb '#555555', \
  "${infile}" using 1:(column(5)/dscale) title "(${kind})'/${dscale}" with lines lw 0.5 lc rgb '#ff8800', \
  "${infile}" using 1:(column(6)/dscale) title "(generic)'/${dscale}" with lines lw 0.5 lc rgb '#008888', \
  "${infile}" using 1:(0.5) notitle with lines lt 0
quit
EOF

atril ${epsfile}
