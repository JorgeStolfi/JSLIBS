#! /bin/bash
# Last edited on 2023-12-14 11:11:11 by stolfi

# Tests the proper weights for unidimensional Shepard-like interpolation of first order.

tmp="/tmp/$$"

color=( 'ff0000' 'aa3300' '557700' '008844' '0033ff' '7700ff' )
show="SHOW"

ctrfile="out/centers.txt"
wtfile="out/weights.txt"
basfile="out/basis.txt"
pngfile="out/plot.png"

generate_files.gawk

# Get the number of centers:
nc=`cat ${ctrfile} | wc -l`

# Create the plot command files:
sep=""
i=0 # Index of center point.
icolor=0 # Index of line color.

# Weight plotting command:
wtgpl="${tmp}-w.gpl"
printf "plot '${ctrfile}' using 2:(0.0) title 'centers' with points lt 1 lw 1.5 pt 7 ps 2.0 lc rgb '#0077ff'" > ${wtgpl}

# Basis plotting command:
basgpl="${tmp}-b.gpl"
printf "plot '${ctrfile}' using 1:(0.0) title 'centers' with points lt 1 lw 3.0 pt 7 ps 2.0 lc rgb '#0077ff', "'\\'"\n" > ${basgpl}
printf "     ''           using 1:(1.0) title 'data'    with points lt 1 lw 3.0 pt 7 ps 2.0 lc rgb '#0077ff'" >> ${basgpl}

while [[ ${i} -lt ${nc} ]]; do
  printf ", "'\\'"\n" >> ${wtgpl}
  printf "     '${wtfile}' using 1:(column(2+${i})) title 'wt[${i}]' with lines lt 1 lw 3.0 lc rgb '#${color[${icolor}]}'" >> ${wtgpl}

  printf ", "'\\'"\n" >> ${basgpl}
  printf "     '${basfile}' using 1:(column(2+${i})) title 'bas[${i}]' with lines lt 1 lw 3.0 lc rgb '#${color[${icolor}]}'" >> ${basgpl}
  
  i=$(( ${i} + 1 ));
  icolor=$(( ${icolor} + 1 ));
  if [[ ${icolor} -ge ${#color[@]} ]]; then icolor=0; fi
done

printf "\n" >> ${wtgpl}

printf "\n" >> ${basgpl}

export GDFONTPATH=.:${HOME}/ttf
gnuplot << EOF
set term png size 2800,2000 font "arial,20"
set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 3.0 lc rgb '#ffddaa'
set output "${tmp}.png"
set multiplot layout 1,1 title "First-order Shepard interpolation"
yden = 2000.0/1800.0
# ----------------------------------------------------------------------
# Pairs
set origin 0.0,(0.500/yden)
set size 1.0,(0.475/yden)
set yrange [-0.01:3.01]
load "${wtgpl}"
# ----------------------------------------------------------------------
# Difference
set origin 0.0,(0.000/yden)
set size 1.0,(0.475/yden)
set yrange [-3.01:3.01]
load "${basgpl}"
# ----------------------------------------------------------------------
unset multiplot
quit
EOF

convert ${tmp}.png -resize '50%' ${pngfile}

if [[ "/${show}" == "/SHOW" ]]; then
  display ${pngfile}
fi

rm -fv ${wtgpl} ${basgpl} ${tmp}.*
