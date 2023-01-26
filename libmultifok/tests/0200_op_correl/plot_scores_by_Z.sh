#! /bin/bash
# Last edited on 2023-01-25 07:56:59 by stolfi

# Reads the file {cdataFile} = "{inPrefix}-cdata.txt" which is supposed to contain lines
# 
#   "P{ki}.{ix}.{iy} {wp} {sharp} {score} {zave} {zdev}"
#
# where {ki} is the input image index, {ix,iy} are column and row of the
# pixel, {wp} is a pixel weight, {sharp} is the mean sharpness in the
# pixel, {score} is the sharpness score computed by
# {multifok_focus_op_score_from_basis}, {zave} is the average scene {Z}
# in the pixel relative to the focus plane {Z}, and {zdev} is the
# deviation of the {Z} coord in the pixel.

prefix="$1"; shift      # File name minus the "-cdata.txt" tail.

tmp="/tmp/$$"

cdataFile="${prefix}-cdata.txt"  # File with compute score and other data.

tdataFile="${tmp}-ndata.txt" # File with processed data for plot.

# Extract the {sharp}, {score}, and {zave}, normalize the scores for each image so that 
# the max score is 1.0, write to temp file with blank lines as separators:

wc -l ${cdataFile} 1>&2

cat ${cdataFile} \
  | process_scores_for_plot.gawk \
  > ${tdataFile}

export GDFONTPATH="${HOME}/tt-fonts"

gnuplot << EOF

set term X11 size 600,800

set xlabel "avg rel Z in pixel"

# Fields from the ${cdataFile}:
t_zave(dum) = column(1)
t_sharp(dum) = column(2)
t_score(dum) = column(3)
t_err(dum) = column(4)

set ylabel "norm avg score"
set title "Normalized scores by relative Z"

plot \
  "${tdataFile}" using (t_zave(0)):(t_score(0)) title "normalized score" with linespoints pt 7 ps 0.75 lc rgb '#ff0000'

pause mouse button1

set ylabel "score - sharp"
set title "Score error by relative Z"

plot \
  "${tdataFile}" using (t_zave(0)):(t_err(0)) title "score-sharp" with linespoints pt 7 ps 0.75 lc rgb '#ff0000'

pause mouse button1

EOF
