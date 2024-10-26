#! /bin/bash 
# Last edited on 2024-10-26 07:59:39 by stolfi

# Plots the various image values of a selected pixel as a function of {zFoc}. 
#
# Reads a file "{dataFile}.txt" produced by
# {test_mfok_scene_make_frame.c} with one line per frame image,
# containing
#
#   {ki} {zFoc} {zDep}  {hAvg} {hDev} {shrp} {sVal[0]} ...  {sVal[NC-1]}
#
# where {ki} is a frame index in the stack, {zFoc} is the {Z}-coordinate
# of the in-focus plane, {zDep} is the nominal depth of focus, {hAvg[k]}
# and {hDev[k]} are the average and deviation of the scene's height
# within that pixel, {shrp[k]} is the pixel's sharpness indicator
# (1/{vBlr}, where {vBlr} is the average of the squared horizontal
# radius of the hit points), and {sVal[0..NC-1]} are the values of
# channels {0..NC-1} of the simulated scene view {sVal}.
#
# Creates an image "{dataFile}.png" with plots of {hAvg} {hDev} {shrp} {sVal[0..NC-1]}
# as a function of {zFoc}.

dataFile="$1"; shift # File with pixel data, with the ".txt" extension.

jobName="${dataFile/.txt/}"

# Determine number ${ni} of frames in stack:
ni=`cat ${dataFile} | egrep -e '^[ ]*[-+0-9.]' | wc -l`
echo "found ${ni} frames in data file" 1>&2 

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
      pcols="(zrav(${k})):(hrad2(${k})):(hDev(${k})) notitle"
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
set xlabel "zFoc"
set ylabel "value"
set xtics 5.0
set mxtics 5
set ytics 0.1
set mytics 5
set yrange [0:3.0]

set grid xtics ytics mxtics lt 3 lw 1.5 lc rgb '#ffddaa', lt 0 lw 1 lc rgb '#ffddaa'

# Scene coordinate range: */
zMin = 0.0
zMax = 30.0

znorm(z) = (z - zMin)/(zMax - zMin)

# {ki} {zFoc} {zDep} {hAvg} {hDev} {shrp} {sVal.R} {sVal.G} {sVal.B}
zFoc(k) = column(2)
zDep(k) = column(3)
hAvg(k) = znorm(column(4))
hDev(k) = column(5)/(zMax - zMin)
shrp(k) = column(6)
sRed(k) = column(7)
sGrn(k) = column(8)
sBlu(k) = column(9)

zUnc(k) = zDep(k)
hDif(k) = column(4) - zFoc(k)
vBlr(k) = 0.5/shrp(k)

plot \
  "${dataFile}" using (zFoc(0)):(shrp(0)):(zUnc(0)) title "zDep" with xerrorbars lc rgb '#777777', \
  ""            using (zFoc(0)):(hAvg(0)) title "hAvg" with linespoints pt 7 ps 2.0 lc rgb '#883344', \
  ""            using (zFoc(0)):(hDev(0)) title "hDev" with linespoints pt 7 ps 2.0 lc rgb '#cc7788', \
  ""            using (zFoc(0)):(shrp(0)) title "shrp" with linespoints pt 7 ps 2.0 lc rgb '#335577', \
  ""            using (zFoc(0)):(vBlr(0)) title "vBlr" with linespoints pt 7 ps 2.0 lc rgb '#995522', \
  ""            using (zFoc(0)):(sRed(0)) title "sVal.R" with linespoints pt 7 ps 2.0 lc rgb '#ff2200', \
  ""            using (zFoc(0)):(sGrn(0)) title "sVal.G" with linespoints pt 7 ps 2.0 lc rgb '#008800', \
  ""            using (zFoc(0)):(sBlu(0)) title "sVal.B" with linespoints pt 7 ps 2.0 lc rgb '#0055ff'

pause mouse button1
EOF

if [[ -s ${tempImage} ]]; then
  convert ${tempImage} -resize '50%' ${plotImage}
  display -title '%f' ${plotImage}
else
  echo "** plot failed" 1>&2 ; exit 1
fi

rm ${tmp}-*
