#! /bin/bash
# Last edited on 2025-03-16 17:43:59 by stolfi

ifile="fourob-MF-0512x0422-Z.fni"
pfile="fourob-MF-0512x0422-Z-scatter.png"

tmp="/tmp/$$"
tfile="${tmp}.png"

export GDFONTPATH="${HOME}/ttf"
gnuplot <<EOF

set term png size 1600,1200 font "arial,24"
set output "${tfile}"

set xlabel "height Z"

plot "${ifile}" using 3:4 title "weight" with points pt 7 ps 1.0 lc rgb '#ff00ff'

EOF

if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile}
  display ${pfile}
fi

hfile_0="${ifile/.fni/}-0-hist.eps"
fni_hist -step 0.02 -channel 0 -logScale < ${ifile} > ${hfile_0} ; evince ${hfile_0}

hfile_1="${ifile/.fni/}-0-hist.eps"
fni_hist -step 0.02 -channel 1 -logScale < ${ifile} > ${hfile_1} ; evince ${hfile_1}
