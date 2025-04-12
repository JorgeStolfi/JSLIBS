#! /bin/bash 
# Last edited on 2025-03-15 23:02:38 by stolfi

tmp="/tmp/$$"

test_interpolate -outDir out

export GFONTPATH=ttf
 
tfile="${tmp}.png" 
gnuplot <<EOF
set term png size 3200,1600 font "arialb,18"
set output "${tfile}"

set xrange [-0.7:]

set multiplot
# ----------------------------------------------------------------------  
set origin 0.000,0.350; set size 1.000,0.650
set lmargin 10
plot \
  "out/data.txt"   using 1:2 title "data" with points pt 7 ps 1.5 lc rgb '#0044ff', \
  ""               using 1:3 notitle with impulses lw 2 lc rgb '#008833', \
  ""               using 1:3 notitle with points pt 7 ps 1.0 lc rgb '#008833', \
  "out/interp.txt" using 1:2 title "expect" with lines lc rgb '#4466ff', \
  ""               using 1:3 title "interp" with points pt 7 ps 1.0 lc rgb '#ff2200', \
  ""               using 1:4 title "weight" with impulses lw 2 lc rgb '#995500', \
  ""               using 1:4 notitle with points pt 7 ps 1.0 lc rgb '#995500'
# ----------------------------------------------------------------------  
set origin 0.000,0.000; set size 1.000,0.350
set lmargin 10

szlog(zz) = log(sqrt(zz*zz + 1) + zz + 1) - log(sqrt(zz*zz + 1) - zz +1)
slog(tt) = szlog(tt/1.0e-9)
slogk(k) = slog(column(k))

set ytics ( \
       "-1"    (slog(-1.0)) 0, \
    "-1e-1" (slog(-1.0e-1)) 0, \
    "-1e-2" (slog(-1.0e-2)) 0, \
    "-1e-3" (slog(-1.0e-3)) 0, \
    "-1e-4" (slog(-1.0e-4)) 0, \
    "-1e-5" (slog(-1.0e-5)) 0, \
    "-1e-6" (slog(-1.0e-6)) 0, \
    "-1e-7" (slog(-1.0e-7)) 0, \
    "-1e-8" (slog(-1.0e-8)) 0, \
    "-1e-9" (slog(-1.0e-9)) 1, \
        "0"       (slog(0)) 0, \
    "+1e-9" (slog(+1.0e-9)) 1, \
    "+1e-8" (slog(+1.0e-8)) 0, \
    "+1e-7" (slog(+1.0e-7)) 0, \
    "+1e-6" (slog(+1.0e-6)) 0, \
    "+1e-5" (slog(+1.0e-5)) 0, \
    "+1e-4" (slog(+1.0e-4)) 0, \
    "+1e-3" (slog(+1.0e-3)) 0, \
    "+1e-2" (slog(+1.0e-2)) 0, \
    "+1e-1" (slog(+1.0e-1)) 0, \
       "+1"    (slog(+1.0)) 0 \
  )
set xzeroaxis 
set yrange [ slog(-1.1e0) : slog(+1.1e0) ]

plot \
  "out/interp.txt" using 1:(slogk(5)) notitle       with impulses lw 2 lc rgb '#aa2200', \
  ""               using 1:(slogk(5)) title "error" with points pt 7 ps 1.0 lc rgb '#aa2200'

# # ----------------------------------------------------------------------
# set origin 0.100,0.300; set size 0.400,0.400
# 
# set xrange [-1.0e-5 : +1.0e-5]
# plot slog(x) with lines

unset multiplot
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

