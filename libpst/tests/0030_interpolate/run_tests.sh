#! /bin/bash 
# Last edited on 2025-02-20 06:55:15 by stolfi

tmp="/tmp/$$"

test_interpolate -outDir out

export GFONTPATH=ttf
 
tfile="${tmp}.png" 
gnuplot <<EOF
set term png size 3200,1200 font "arialb,18"
set output "${tfile}"

plot \
  "out/data.txt"   using 1:2 title "data" with points pt 7 ps 1.0 lc rgb '#0044ff', \
  ""               using 1:3 notitle with impulses lw 2 lc rgb '#008833', \
  ""               using 1:3 notitle with points pt 7 ps 1.0 lc rgb '#008833', \
  "out/interp.txt" using 1:2 title "interp" with points pt 7 ps 1.0 lc rgb '#ff2200', \
  ""               using 1:3 notitle with impulses lw 2 lc rgb '#995500', \
  ""               using 1:3 notitle with points pt 7 ps 1.0 lc rgb '#995500'
EOF

if [[ -s ${tfile} ]]; then
  pfile="out/plot.png"
  convert ${tfile} -resize '50%' ${pfile}
  display ${pfile}
else
  echo "** file ${tfile} not generated" 1>&2; exit 1
fi

rm -f ${tfile}

echo "OK" 1>&2

