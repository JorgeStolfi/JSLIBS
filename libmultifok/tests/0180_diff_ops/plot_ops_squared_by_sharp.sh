#! /bin/bash
# Last edited on 2023-01-25 18:29:57 by stolfi

# Reads the file {odataFile} = "{prefix}-odata.txt" which is supposed to contain lines
# 
#   "P{ki}.{ix}.{iy} {vave} {vdev} {wp} {sharp} {zave} {zdev} {coeff[0]} .. {coeff[NB-1]}"
#
# where {ki} is the input image index, {ix,iy} are column and row of the
# pixel, {wp} is a pixel weight, {sharp} is the mean sharpness in the
# pixel, {score} is the sharpness score computed by
# {multifok_focus_op_score_from_basis}, {zave} is the average scene {Z}
# in the pixel relative to the focus plane {Z}, and {zdev} is the
# deviation of the {Z} coord in the pixel.
#
# Writes a file ${chistFile} = "${prefix}-chist.txt" with the avreage coeffs squared per 
# sharpness bin. See {process_coeffs_for_plot.gawk}.

prefix="$1"; shift      # File name minus the "-odata.txt" tail.
basisName="$1"; shift   # Basis name ("DIFF", "LAPL", "HART", etc).
unitTerm="$1"; shift    # True (1) to include the unit term in the regression.
title="$1"; shift       # Plot title.

# Input files:
odataFile="${prefix}-odata.txt"       # File with compute score and other data.
belNameFile="${prefix}-belnames.txt"  # File with basis element names.

# Output files:
plotImageFile="${prefix}-cfhist.png"             # Plot of binned coeffs squared by {sharp}. 
chistFile="${prefix}-un${unitTerm}-chist.txt"    # File with histogram data for plot.
outFormFile="${prefix}-un${unitTerm}-oform.txt"  # Fitted formula for sharp as function of the coeffs squared.
outRegrFile="${prefix}-un${unitTerm}-oregr.txt"  # File with true and fitted sharpness.
outCombFile="${prefix}-un${unitTerm}-ocform.txt" # Fitted formula with related terms equalized.

# Get the basis element names and count:
belName=( `cat ${belNameFile}` )
nb=${#belName[@]}
if [[ ${nb} -gt 9 ]]; then echo "** too many coeffs" 1>&2; exit 1; fi
echo "found ${nb} basis elements" 1>&2 
termName=( `cat ${belNameFile} | sed -e 's:^\([A-Z][A-Z0-9]*\)$:\1*\1:g' ` )
echo "termName = ${termName[*]}" 1>&2 

tmp="/tmp/$$"

tmpGplFile="${tmp}-plot.gpl"    # File with gnuplot commands.
tmpImageFile="${tmp}-plot.png" # Temporary plot file.

wc -l ${odataFile} 1>&2

cat ${odataFile} \
  | ./process_coeffs_for_plot.gawk \
      -v nb=${nb} \
      -v unitTerm=${unitTerm} \
  > ${chistFile}
wc -l ${chistFile}
  
# Assemble the plot command:
color=( '#cc7700' '#888800' '#008800' '#007755' '#0033aa' '#0000ff' '#5500ff' '#aa0088' '#555555' '#ff0000' )
kb=0
printf "plot" > ${tmpGplFile}
sep=""
while [[ ${kb} -lt ${nb} ]]; do
  printf "  element bas[%d] = %s...\n" ${kb} "${termName[${kb}]}"
  printf " "${sep}'\\'"\n" >> ${tmpGplFile}
  printf "  \"${chistFile}\" using 2:(column(4+${kb})) title \"${termName[${kb}]}\"" >> ${tmpGplFile}
  printf " with linespoints lw 2 pt 7 ps 1.5 lc rgb '${color[${kb}]}'"  >> ${tmpGplFile}
  sep=","
  kb=$(( ${kb} + 1 ))
done
printf "\n" >> ${tmpGplFile}

echo "plot commands:" 1>&2 
cat  ${tmpGplFile} 1>&2
  
export GDFONTPATH="${HOME}/tt-fonts"
rm -f ${plotImageFile}
gnuplot << EOF
set term png size 1600,1500 noenhanced font "arial,24"
set output "${tmpImageFile}"

set title "${title}"
set xlabel "pixel sharpness"
set ylabel "mean squared basis coefficient"

load "${tmpGplFile}"

pause mouse button1

EOF

if [[ -s ${tmpImageFile} ]]; then
  convert ${tmpImageFile} -resize '50%' ${plotImageFile}
  display -title '%f' ${plotImageFile}
else
  echo "** plot failed" 1>&2 ; exit 1
fi

# Do a regression on the histogram file:

if [[ ${unitTerm} -ne 0 ]]; then
  termName+=( "1" )
fi

nt=${#termName[@]}

linear_fit \
    -terms ${nt} \
    -weighted T \
    -termNames "${termName[@]}" \
    -writeFormula ${outFormFile} \
  < ${chistFile} \
  > ${outRegrFile}

plot_regression_result.sh SHOW ${outRegrFile/.txt/}
combine_${basisName}_coeffs.gawk < ${outFormFile} > ${outCombFile}; cat ${outCombFile} 1>&2
 
rm ${tmp}-*
