#! /bin/bash 
# Last edited on 2023-04-26 14:55:55 by stolfi

# Plots the sharpness as a function of Z for selected pixels 
# in an frame stack.
#
# Reads a file "{fileName}.txt" produced by {mf_0100_make_image.c} with one line per frame image,
# containing
#
#   {zFoc} {zDep} {zave[0]} {zdev[0]} {sharp[0]} ...  {zave[nq-1]} {zdev[nq-1]} {sharp[nq-1]}
#
# where {zFoc} is the {Z}-coordinate of the focused plane, {zDep} is the depth of focus,
# and each triplet {zave[k]} {zdev[k]} {sharp[k]} is the average and deviation of the 
# scene {Z} coordinate in some selected pixel, and {sharp[k]} is the pixel's sharpness (the mean radius of 
# the sampling kernel). 
#
# Creates an image "{fileName}.png" with {nq} plots, for each selected pixel, with {sharp[k]} 
# as a function of {zave[k]-zFoc}.

dataFile="$1"; shift # File with pixel data, with the ".txt" extension.

jobName="${dataFile/.txt/}"

# Determine the number of pixels ${nq}:
fld=( `cat ${dataFile} | egrep -e '^[ ]*[-+0-9.]' | head -n 1` )
nf=${#fld[@]}
if [[ ${nf} -eq 2 ]]; then echo "!! no pixels in file" 1>&2 ; exit 0; fi
if [[ ${nf} -le 2 ]]; then echo "** bad number of fields - ${nf} fields" 1>&2 ; exit 1; fi
nq=$(( ( ${nf} - 2) / 3 )) # Number of pixels in file.
if [[ ${nf} -ne $(( 2 + 3 * ${nq} )) ]]; then echo "** bad number of fields = ${nf}" 1>&2 ; exit 1; fi

# Determine number ${ni} of frames in stack:
ni=`cat ${dataFile} | egrep -e '^[ ]*[-+0-9.]' | wc -l`
echo "found ${ni} frames and ${nq} pixels in file" 1>&2 

tmp="/tmp/$$"

gplFile="${tmp}.gpl"

# Assemble the plot command:
color=( '#ff0000' '#cc7700' '#888800' '#008800' '#007755' '#0033aa' '#9999ff' '#dd77ff' '#7700ff' '#aa00cc' '#555555' )
printf "plot" > ${gplFile}
sep=""
for style in xerrorbars linespoints ; do 
  k=0
  kc=0
  while [[ ${k} -lt ${nq} ]]; do
    if [[ "/${style}" == "/xerrorbars" ]]; then
      pcols="(zrav(${k})):(hrad2(${k})):(zdev(${k})) notitle"
      szops="pt 0 lw 3"
    else
      pcols=`printf "(zrav(${k})):(hrad2(${k})) title \"%02d\"" "${k}"`
      szops="pt 7 ps 2.0 lw 3"
    fi
    printf "${sep} "'\\'"\n" >> ${gplFile}
    printf "  \"${dataFile}\" using ${pcols}" >> ${gplFile}
    printf " with ${style} ${szops} lc rgb '${color[kc]}'" >> ${gplFile}
    sep=","
    k=$(( ${k} + 1 ))
    kc=$(( ${kc} + 1 ))
    if [[ ${kc} -ge ${#color[@]} ]]; then kc=0; fi
  done
done
printf "\n"

tempImage="${tmp}-d.png"
plotImage="${jobName}.png"

export GDFONTPATH="${HOME}/tt-fonts"

gnuplot << EOF
set term png size 2800,1500 noenhanced font "arial,20"
set output "${tempImage}"

set key top left reverse Left
set xlabel "zave[k] - zFoc"
set ylabel "hrad2[k]"
set xtics 5.0
set mxtics 5
set ytics 5
set mytics 5

set grid xtics ytics mxtics lt 3 lw 1.5 lc rgb '#ffddaa', lt 0 lw 1 lc rgb '#ffddaa'

zFoc(k) = column(1)
zDep(k) = column(2)
zave(k)= column(3 + 3 * k)
zrav(k)= zave(k) - zFoc(k)
zdev(k)= column(3 + 3 * k + 1)

sharp(k) = (zdev(k) > 1.5 ? 0/0 : column(3 + 3 * k + 2))
hrad(k) = (zdev(k) > 1.5 ? 0/0 : 1.0/sharp(k))
hrad2(k) = (zdev(k) > 1.5 ? 0/0 : hrad(k)**2)

load "${gplFile}"

pause mouse button1
EOF

if [[ -s ${tempImage} ]]; then
  convert ${tempImage} -resize '50%' ${plotImage}
  display -title '%f' ${plotImage}
else
  echo "** plot failed" 1>&2 ; exit 1
fi

rm ${tmp}-*
