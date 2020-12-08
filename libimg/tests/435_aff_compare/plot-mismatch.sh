#! /bin/bash
# Last edited on 2020-10-13 22:16:20 by jstolfi

prefix="$1"; shift
deform="$1"; shift

dfile="out/${prefix}_${deform}_plot.dat"
pfile="out/${prefix}_${deform}_plot.png"

gnuplot <<EOF
  set term png medium size 1200,1200 font "arial,20"
  set output "out/.tmp.png"
  set title "${prefix/_/ } ${deform}"
  set hidden3d
  sel(k,v,i) = (column(k) == v ? column(i) : 0/0)
  splot \
    "${dfile}" using (sel(1,0,4)):5:6 title "mismatch" with lines lw 3, \
    "${dfile}" using (sel(1,1,4)):5:6 title "optimum" with points ps 2
EOF

convert out/.tmp.png -resize '50%' ${pfile} 
display ${pfile}
