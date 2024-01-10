#! /bin/bash 
# Last edited on 2024-01-05 18:12:34 by stolfi

dfile="$1"; shift; # File with Fourier coeffs before and after cleanup.
Gmax="$1"; shift;  # Mag gain modulus for plot scales.

pfile="${dfile/.txt/}.png"
rm -f ${pfile}
  
tmp="/tmp/$$"
tfile="${tmp}.png"

export GDFONTPATH=.:${HOME}/ttf
gnuplot <<EOF
set term png size 1600,1600 font "arial,24"
set output "${tfile}"

Gmin = 1.0e-16
Gmax = ${Gmax}
set logscale x
set xrange [(Gmin):(Gmax)]
set yrange [(Gmin):(Gmax)]

set logscale y

set grid xtics
set grid ytics

gain(k) = (column(k) == 0 ? 0/0 : column(k))
gratio(i,j) = gain(i)/gain(j)

plot \
  "${dfile}" using (gain(4)):(gain(7)) title "|F'| : |F|" with linespoints lw 2 pt 7 ps 2 lc rgb '#338800'
EOF

if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile}
  rm ${tfile}
  if [[ -s ${pfile} ]]; then
    display -title "${dfile/.txt/}" ${pfile}
  else
    echo "** resize failed" 1>&2; exit 1
  fi
else
  echo "** gnuplot failed" 1>&2; exit 1
fi

