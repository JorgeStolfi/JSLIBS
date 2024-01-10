#! /bin/bash 
# Last edited on 2024-01-04 04:12:25 by stolfi

dfile="$1"; shift; # File with tabulated Fourier and Hartley gains. 
fsmp="$1"; shift;  # Nominal sampling frequency (Hz).

pfile="${dfile/.txt/.png}"
rm -f ${pfile}
  
tmp="/tmp/$$"

# Get min and max freq in data file:
fmin=`head -n 10 ${dfile} | egrep -e '^ *[#] *fmin *[=]' | sed -e 's:^.*[=] *::g'`
imin=`gawk "BEGIN{ imin = int((-(${fmin}))/${fsmp} + 0.999999); print  imin}"`
echo "fmin = ${fmin} imin = ${imin}" 1>&2

fmax=`head -n 10 ${dfile} | egrep -e '^ *[#] *fmax *[=]' | sed -e 's:^.*[=] *::g'`
imax=`gawk "BEGIN{ imax = int((-(${fmax}))/${fsmp} + 0.999999); print  imax}"`
echo "fmax = ${fmax} imax = ${imax}" 1>&2

# Create file with key frequencies:
ifile="${tmp}_i.txt"
rm -f ${ifile}
nimp=10
i=${imin}
while [[ ${i} -le ${imax} ]]; do
  fsi=`echo "${i}*${fsmp}" | bc -lq`  # Sampling frequency.
  fni=`echo "(${i} + 0.5)*${fsmp}" | bc -lq`  # Nyquist frequency.
  echo "${fni} -100.0" >> ${ifile}
  echo "${fni} +100.0" >> ${ifile}
  echo "${fsi} -100.0" >> ${ifile}
  echo "${fsi} +100.0" >> ${ifile}
  i=$(( $i + 1 ))
done

tfile="${tmp}.png"
export GDFONTPATH=.:${HOME}/ttf
gnuplot <<EOF

set term png size 2800,1200 font "arial,24"
set output "${tfile}"
fmin = ${fmin}
fmax = ${fmax}
df = fmax-fmin
logy = 0
LGmax = 9
set grid xtics
set grid ytics
set xrange [(fmin-0.05*df):(fmax+0.05*df)]

signlog(y) = (y == 0 ? 0 : (y < 0 ? -1 : +1)*hidesmall(safelog(abs(y))))
safelog(y) = log(y)/log(10) + LGmax + 1
hidesmall(s) = (s < 1 ? 0/0 : s)

if (logy) \
  { 
    yval(k) = signlog(column(k))
    set yrange [(-1.2*LGmax):(+1.2*LGmax)]
    set ytics ( \
      "-1"      -10 0, \
      "-0.1"     -9 0, \
      "-0.01"    -8 0, \
      "-0.001"   -7 0, \
      "-0.0001"  -6 0, \
      "-1.0e-5"  -5 0, \
      "-1.0e-6"  -4 0, \
      "-1.0e-7"  -3 0, \
      "-1.0e-8"  -2 0, \
      "-1.0e-9"  -1 0, \
      "0"         0 0, \
      "+1.0e+9"   1 0, \
      "+1.0e+8"   2 0, \
      "+1.0e+7"   3 0, \
      "+1.0e+6"   4 0, \
      "+1.0e+5"   5 0, \
      "+0.0001"   6 0, \
      "+0.001"    7 0, \
      "+0.01"     8 0, \
      "+0.1"      9 0, \
      "+1"       10 0  \
    ) 
  } \
else \
  { 
    yval(k) = (column(k) < -50 ? 0/0 : column(k))
    set yrange [-1.6:+1.6]
  }

plot \
  "${ifile}" using 1:2 notitle with impulses lw 6 lc rgb '#ffccaa', \
  "${dfile}" using 2:(yval(3)) title "re(F) org" with linespoints pt 7 ps 1 lw 3 lc rgb '#0022ff', \
  "${dfile}" using 2:(yval(4)) title "im(F) org" with linespoints pt 7 ps 1 lw 3 lc rgb '#ff0033', \
  "${dfile}" using 2:(yval(5)) title "H raw"     with linespoints pt 7 ps 1 lw 3 lc rgb '#008833', \
  "${dfile}" using 2:(yval(6)) title "H kuk"     with linespoints pt 7 ps 1 lw 3 lc rgb '#aa6600', \
  "${dfile}" using 2:(yval(7)) title "re(F) rec" with linespoints pt 6 ps 2 lw 1 lc rgb '#88aaff', \
  "${dfile}" using 2:(yval(8)) title "im(F) rec" with linespoints pt 6 ps 2 lw 1 lc rgb '#ff88aa'
EOF

rm -f ${ifile}
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

