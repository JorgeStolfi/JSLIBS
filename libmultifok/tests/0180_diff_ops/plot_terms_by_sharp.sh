#! /bin/bash
# Last edited on 2023-02-01 18:24:24 by stolfi

# Reads the file {histDataFile} = "{prefix}-rhdata.txt" which is supposed to contain one line
# for each sharpness bin, with data
# 
#   "{kh} {sharp} {term[0..nt-1]}"
#
# where {kk} is the sharpness bin index, {sharp} is the central sharpness of that bin, in the
# pixel, and {term[0..nt-1]} are the average values of the {nt} quadratic terms to use 
# for regression.

echo "=== plot_terms_by_sharp.sh =============================" 1>&2
echo "$@" 1>&2

prefix="$1"; shift      # File name minus the "-odata.txt" tail.
basisName="$1"; shift   # Basis name ("DIFF", "LAPL", "HART", etc).
termSet="$1"; shift     # Term set ("ALL", etc.).
title="$1"; shift       # Plot title.

# Input files:
belNameFile="${prefix}-bnames.txt"      # File with basis element names.
termNameFile="${prefix}-tnames.txt"     # File with quadratic term names.

# Output files:
histDataFile="${prefix}-hdata.txt"      # File with histogram data for plot and regression.
histPlotFile="${prefix}-hdata.png"      # Plot of binned coeffs squared by {sharp}. 

# Get the basis element names and count:
belName=( `cat ${belNameFile}` )
nb=${#belName[@]}
echo "found ${nb} basis elements" 1>&2 
if [[ ${nb} -gt 9 ]]; then echo "** too many coeffs" 1>&2; exit 1; fi

# Get the term names and count:
termName=( `cat ${termNameFile}` )
nt=${#termName[@]}
echo "found ${nt} quadratic terms" 1>&2 
if [[ ${nt} -gt $(( ${nb} * (${nb} + 1) / 2 )) ]]; then echo "** too many terms" 1>&2; exit 1; fi

tmp="/tmp/$$"

tmpGplFile="${tmp}-plot.gpl"    # File with gnuplot commands.
tmpPlotFile="${tmp}-plot.png" # Temporary plot file.
  
# Assemble the plot command:
color=( '#ff0000' '#cc7700' '#888800' '#008800' '#007755' '#0033aa' '#9999ff' '#dd77ff' '#7700ff' '#aa00cc' '#555555' )
kt=0
printf "plot" > ${tmpGplFile}
sep=""
while [[ ${kt} -lt ${nt} ]]; do
  printf "  term[%d] = %s...\n" ${kt} "${termName[${kt}]}"
  printf " "${sep}'\\'"\n" >> ${tmpGplFile}
  tnk="${termName[${kt}]}"
  printf "  \"${histDataFile}\" using 2:(column(3+${kt})) title \"term[%02d]\"" ${kt} >> ${tmpGplFile}
  printf " with linespoints lw 2 pt 7 ps 2.0 lc rgb '${color[${kt}]}'"  >> ${tmpGplFile}
  sep=","
  kt=$(( ${kt} + 1 ))
done
printf "\n" >> ${tmpGplFile}

echo "plot commands:" 1>&2 
cat  ${tmpGplFile} 1>&2
  
export GDFONTPATH="${HOME}/tt-fonts"
rm -f ${histPlotFile}
gnuplot << EOF
set term png size 1600,1500 noenhanced font "arial,24"
set output "${tmpPlotFile}"

set title "${title}"
set xlabel "pixel sharpness"
set ylabel "quadratic term value"

set key top left

load "${tmpGplFile}"

pause mouse button1

EOF

if [[ -s ${tmpPlotFile} ]]; then
  convert ${tmpPlotFile} -resize '50%' ${histPlotFile}
  display -title '%f' ${histPlotFile}
else
  echo "** plot failed" 1>&2 ; exit 1
fi

rm ${tmp}-*
