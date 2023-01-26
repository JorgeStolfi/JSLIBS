#! /bin/bash
# Last edited on 2023-01-25 07:55:22 by stolfi

# Reads file {regrFile} = "{inPrefix}-regr.txt" which is supposed to contain lines
# 
#   "P{ki}.{ix}.{iy} {sharp^exp} {sregr^exp}"
#
# where {ki} is the input image index, {ix,iy} are column and row of the
# pixel, where {sharp} is the mean sharpness in
# the pixel, {sregr} is the result of the regression fitted fotmula,
# and {exp} is 2 if {sqrSharp} is true else 1.

# Also reads the file {cdataFile} = "{inPrefix}-cdata.txt" which is supposed to contain lines
# 
#   "P{ki}.{ix}.{iy} {wp} {sharp} {score} {zave} {zdev}"
#
# where {wp} is a pixel importance weight, {score} is the sharpness
# score computed by {multifok_focus_op_score_from_basis}, {zave} is the
# average scene {Z} in the pixel, relative to the focus plane {Z}, and
# {zdev} is the deviation of the {Z} coord in the pixel.

prefix="$1"; shift      # File name minus the "-regr.txt" or "-cdata.txt" tail.
sqrSharp="$1"; shift    # If 1, assume that {sharp} and {score} are squared in {regrFile}.
title="$1"; shift       # Plot title.

tmp="/tmp/$$"

regrFile="${prefix}-regr.txt"  # Regression result file.
cdataFile="${prefix}-cdata.txt"  # File with actual and computed sharpness and term coefs.

if [[ ${sqrSharp} -ne 0 ]]; then xsqr=" squared"; else xsqr=""; fi

export GDFONTPATH="${HOME}/tt-fonts"

# Join the files by pixel index:
rfile="${tmp}-regr.txt"
cat ${regrFile} | sort -b -k1 > ${rfile}

cfile="${tmp}-cdata.txt" 
cat ${cdataFile} | sort -b -k1 > ${cfile}

# Fields: "{pixelID} {weight} {comp_sharp} {comp_score} {regr_sharp^exp} {regr_score^exp}"
jfile="${tmp}-join.txt" 
join -j 1 -a 1 -a 2 -e '??' -o0,2.2,2.3,2.4,1.2,1.3 ${rfile} ${cfile} > ${jfile}

# Combine the data into bins:
regrPlotFile="${tmp}-regr.txt"
compPlotFile="${tmp}-true.txt"  

function bin_data() {
  col="$1"; shift;   # Column to bin
  ofile="$1"; shift; # Output file 
  cat ${jfile} \
    | compute_score_av_by_sharp.gawk \
        -v sqrSharp=${sqrSharp} \
        -v col=${col} \
        -v nh=50 \
    > ${ofile}
}

bin_data 4 ${compPlotFile}
bin_data 6 ${regrPlotFile}

gnuplot << EOF

set term X11 size 800,800
set title "${title}"

set xrange [-0.01:+2.01]
set xlabel "actual sharpness"

set yrange [-0.01:+2.01]
set ylabel "fitted or computed score"

sqrsh = ${sqrSharp}

# Fields from the ${cdataFile}:
val(k) = column(k)

# Fields from the ${regrPlotFile}:
err(k1,k2) = val(k1) - val(k2)

plot \
  "${compPlotFile}" using (val(2)):(val(3)) title "comp" with points pt 6 ps 1.00 lc rgb '#ff0000', \
  "${regrPlotFile}" using (val(2)):(val(3)) title "regr" with points pt 7 ps 1.00 lc rgb '#009900'

pause mouse button1

set yrange [-2.01:+2.01]
set ylabel "error of fitted or computed score"

plot \
  "${compPlotFile}" using (val(2)):(err(3,2)) title "comp" with points pt 6 ps 1.00 lc rgb '#ff0000', \
  "${regrPlotFile}" using (val(2)):(err(3,2)) title "regr" with points pt 7 ps 1.00 lc rgb '#009900'

pause mouse button1

EOF
