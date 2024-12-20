#! /bin/bash 
# Last edited on 2024-12-19 17:17:31 by stolfi

PROG_NAME=${0##*/}
PROG_DESC="displays a tabulated light response function."
PROG_HELP=(
  "${PROG_NAME} {FILE_NAME} {FUNC_NAME}"
)

# Displays a tabulated sample mapping function, read from a file
# called "{FILE_NAME}.txt".  Assumes that each line has four numbers
# 
#    {IRAW} {VRAW} {VLIN} {VENC} {IENC} {ZERR}  {MRAW} {MLIN} {MENC}
# 
# where 
#
#    {IRAW} (1) is an integer encoded sample value, in the range {-MAXVAL..MAXVAL}
#              
#    {VRAW} (2) is {IRAW/MAXVAL}, in {[-1_+1]}
#              
#    {VLIN} (3) is the decoded linear sample value {D(vRaw)}.
#              
#    {VENC} (4) is the re-encoded value {E(vLin)}.
#              
#    {ZERR} (5) is the difference {VENC-VRAW} times {MAXVAL}.
#              
#    {IENC} (6) is {VENC} times {MAXVAL} rounded to the nearest integer.
#              
#  {MRAW,MLIN,MENC} (7-9) are {VRAW,VLIN,VENC} in a "signed log" scale 
# 
# Assumes that {VRAW} and {VLIN} range in {[-1 _ +1]}, while {MRAW} and
# {MLIN} range in {[-3 _ +3]} or so.

if [[  $# -lt 2 ]]; then
  echo "not enough arguments" 1>&2;
  echo "usage: ${PROG_HELP[@]}" 1>&2; exit 1;
fi

FILE_NAME="$1"; shift; 
FUNC_NAME="$1"; shift; 

if [[  $# -ne 0 ]]; then
  sep="$2"; shift; shift;
  echo "spurious arguments "'"'"$1"'"...' 1>&2;
  echo "usage: ${PROG_HELP[@]}" 1>&2; exit 1;
fi

FILE_NAME="${FILE_NAME/.txt/}" # Just in case.

export GDFONTPATH=ttf

gnuplot <<EOF
set terminal postscript eps color "TimesRoman" 30
set output "${FILE_NAME}-lin.eps"
set size 1,1.45
set size ratio -1
set title "${FUNC_NAME} (lin)"
set xrange[-1.05:+1.05]
set yrange[-1.05:+1.05]
set key top left
plot \
  "${FILE_NAME}.txt" using 2:3 title '(x,D(x))' with linespoints pt 7 ps 0.30 lc rgb '#ff0000', \
  ""                 using 4:3 title '(E(y),y)' with linespoints pt 7 ps 0.30 lc rgb '#0055ff'
quit
EOF

gnuplot <<EOF
set terminal postscript eps color "TimesRoman" 30
set output "${FILE_NAME}-log.eps"
set size 1,1.45
set size ratio -1
set title "${FUNC_NAME} (log)"
set key top left
plot \
  "${FILE_NAME}.txt" using 7:8 title '(x,D(x))' with linespoints pt 7 ps 0.30 lc rgb '#ff0000', \
  ""                 using 9:8 title '(E(y),y)' with linespoints pt 7 ps 0.30 lc rgb '#0055ff'
quit
EOF

gnuplot <<EOF
set terminal postscript eps color "TimesRoman" 30
set output "${FILE_NAME}-err.eps"
set size 3,1
set title "${FUNC_NAME} error (x - E(D(x)))*M"
set nokey
set logscale x
set logscale y
set xrange[:]
set yrange[:]

# To avoid zero in log plot:
fudge(x) = (x < 0 ? 0/0 : x + 1.0e-10)
egduf(x) = (x > 0 ? 0/0 : -x + 1.0e-10)
pos(x) = (x <= 0 ? 0/0 : x)
plot "${FILE_NAME}.txt" using (pos(column(2))):(fudge(column(5))) title 'pos' with linespoints pt 7 ps 0.30 lc rgb '#0055ff', \
     ""                 using (pos(column(2))):(egduf(column(5))) title 'neg' with linespoints pt 7 ps 0.30 lc rgb '#ff0000'
quit
EOF

atril "${FILE_NAME}-lin.eps" &
atril "${FILE_NAME}-log.eps" &
atril "${FILE_NAME}-err.eps"
