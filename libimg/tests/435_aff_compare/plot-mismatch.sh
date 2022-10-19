#! /bin/bash
# Last edited on 2022-10-19 18:55:29 by stolfi

prefix="$1"; shift
deform="$1"; shift

tmp="/tmp/$$"

tfile="${tmp}.png"
dfile="out/${prefix}_${deform}_plot.dat"
pfile="out/${prefix}_${deform}_plot.png"

title="${prefix//_/ } ${deform}"
echo "title = ${title}" 1>&2

gnuplot <<EOF
  set term png medium size 1200,1200 font "arial,20"
  set output "${tfile}"
  set title "${title}"
  set hidden3d
  splot \
    "${dfile}" using 4:5:6 title "mismatch" with lines lw 3
EOF

if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile} 
  display ${pfile}
  rm -f ${tfile}
else
  echo "** plot failed" 1>&2; exit 1
fi


    
#  sel(k,v,i) = (column(k) == v ? column(i) : 0/0)
#  splot \
#    "${dfile}" using (sel(1,0,4)):5:6 title "mismatch" with lines lw 3, \
#    "${dfile}" using (sel(1,1,4)):5:6 title "optimum" with points ps 2
