#! /bin/bash
# Last edited on 2020-10-16 03:47:25 by jstolfi

prefix="$1"; shift
method="$1"; shift

dfile="${prefix}_${method}.dat"
tfile="/tmp/$$.png"
pfile="${prefix}_${method}.png"

title="${prefix/out\//} ${method}"

if [[ ! ( -s ${dfile} ) ]]; then echo "** file ${dfile} not found" 1>&2 ; exit 1 ; fi

gnuplot <<EOF
  set term png medium size 1200,1200 font "arial,24"
  set output "${tfile}"
  set hidden3d
  splot "${dfile}" using 1:2:3 title "${title}" with lines 
EOF

if [[ ! ( -s ${tfile} ) ]]; then echo "** file ${tfile} not creatd" 1>&2 ; exit 1 ; fi

convert ${tfile} -resize '50%' ${pfile}
display -title '%f' ${pfile}
rm -f ${tfile}
