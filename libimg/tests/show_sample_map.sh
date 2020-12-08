#! /bin/bash 
# Last edited on 2007-11-11 00:58:31 by stolfi

PROG_NAME=${0##*/}
PROG_DESC="displays a tabulated light response function."
PROG_HELP=(
  "${PROG_NAME} {NAME}..."
)

# Displays a tabulated sample mapping function, read from a file
# called "{NAME}.txt".  Assumes that each line has four numbers
# 
#    {VRAW} {VCOR} {VINV} {MRAW} {MCOR} {MINV} 
# 
# where {VRAW} is an input value, {VCOR} is the corresponding
# output value, {VINV} is the inverse-mapped {VCOR}; 
# and {MRAW}, {MCOR}, {MINV} are the same value in an
# "signed logarithmic" scale.
# Assumes that {VRAW} and {VCOR} range in {[-1 _ +1]},
# while {MRAW} and {MCOR} range in {[-3 _ +3]} or so.

name="$1"; shift; 

if [[  $# -ne 0 ]]; then
  sep="$2"; shift; shift;
  echo "spurious arguments "'"'"$1"'"...' 1>&2;
  echo "usage: ${PROG_HELP[@]}" 1>&2; exit 1;
fi

gnuplot <<EOF
set terminal postscript eps color "TimesRoman" 18
set output "${name}-lin.eps"
set size 1,1.45
set size ratio -1
set title "${name} (lin)"
set xrange[-1.05:+1.05]
set yrange[-1.05:+1.05]
set nokey
plot "${name}.txt" using 1:2 with linespoints
quit
EOF

gnuplot <<EOF
set terminal postscript eps color "TimesRoman" 18
set output "${name}-inv.eps"
set size 1,1.45
set size ratio -1
set title "${name} (inv)"
set xrange[-1.05:+1.05]
set yrange[-1.05:+1.05]
set nokey
plot "${name}.txt" using 1:3 with linespoints
quit
EOF

gnuplot <<EOF
set terminal postscript eps color "TimesRoman" 18
set output "${name}-log.eps"
set size 1,1.45
set size ratio -1
set title "${name} (log)"
set xrange[-2.10:+2.10]
set yrange[-2.10:+2.10]
set nokey
plot "${name}.txt" using 4:5 with linespoints
quit
EOF

ghostview "${name}-lin.eps" &
ghostview "${name}-inv.eps" &
ghostview "${name}-log.eps"
