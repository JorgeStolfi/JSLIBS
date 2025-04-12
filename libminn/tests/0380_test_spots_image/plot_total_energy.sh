#! /bin/bash
# Last edited on 2025-03-31 15:03:38 by stolfi

dfile="$1"; shift

dname="${dfile%.*}"

export GDFONTPATH=ttf

tfile="out/.tmp.png"
pfile="${dname}.png"

rm -f ${tfile} ${pfile} 

gnuplot <<EOF
  set term png medium size 1200,1200 font "courbd,18"
  set output "${tfile}"
  set title "total energy"
  set xlabel "d0 (px)"
  set ylabel "d1 (px)"
  set zlabel "energy"
  set hidden3d
  splot "${dfile}" using 1:2:3 with lines 
EOF
 
if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile} 
  display -title '%f' ${pfile} 
else
  echo "** ${pfile} not created" 1>&2; exit 1
fi
