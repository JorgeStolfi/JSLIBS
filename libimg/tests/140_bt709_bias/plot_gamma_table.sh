#! /bin/bash
# Last edited on 2019-04-09 22:03:30 by jstolfi

usage="${0} '{TITLE1}' '{TITLE2}' {TABLE_NAME}"

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
set xrange [-0.05:+1.05]
set yrange [-0.05:+1.05]
set title "${title1}\n${title2}"
set xzeroaxis
set yzeroaxis
set key outside
# test terminal
plot \
  "${infile}" using 1:2 title "BT709" with lines lt 1, \
  "${infile}" using 1:3 title "sc_gm" with lines lt 7, \
  "${infile}" using 1:(0.5 + 20*column(4)) title "diff x 20" with lines lt 3, \
  "${infile}" using 1:(0.2*column(5)) title "(BT709)'/5" with lines lw 0.5 lt 1, \
  "${infile}" using 1:(0.2*column(6)) title "(sc\_gm)'/5" with lines lw 0.5 lt 7, \
  "${infile}" using 1:(0.5) notitle with lines lt 0
quit
EOF

atril ${epsfile}
