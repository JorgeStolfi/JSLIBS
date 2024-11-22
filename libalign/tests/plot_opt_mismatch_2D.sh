#! /bin/bash
# Last edited on 2024-11-08 15:02:33 by stolfi

dfile="$1"; shift

prefix="${dfile/.txt/}"

tfile="/tmp/$$.png"
pfile="${prefix}.png"

title="${dfile/out\/test-/}"

if [[ ! ( -s ${dfile} ) ]]; then echo "** file ${dfile} not found" 1>&2 ; exit 1 ; fi

gnuplot <<EOF
  set term png medium size 1200,1200 font "arial,24"
  set output "${tfile}"
  set title "${title}"
  set hidden3d
  splot "${dfile}" using 1:2:3 title "mismatch" with lines lw 3
EOF

if [[ ! ( -s ${tfile} ) ]]; then echo "** file ${tfile} not creatd" 1>&2 ; exit 1 ; fi

convert ${tfile} -resize '50%' ${pfile}
display -title '%f' ${pfile}
rm -f ${tfile}
