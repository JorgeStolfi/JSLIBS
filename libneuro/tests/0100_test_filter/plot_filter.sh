#! /bin/bash 
# Last edited on 2023-12-18 10:56:34 by stolfi

dfile="$1"; shift; # File with gains per frequency. 
fsmp="$1"; shift;  # Nominal sampling frequency (Hz).
flo0="$1"; shift;  # Nominal lowest pass freq (Hz).
flo1="$1"; shift;  # Nominal lowest preserved frew (Hz).
fhi1="$1"; shift;  # Nominal highest preserved freq (Hz).
fhi0="$1"; shift;  # Nominal highest pass frew (Hz). 

pfile="${dfile/.txt/.png}"
rm -f ${pfile}
  
tmp="/tmp/$$"

# Determine number of curves {nc}:
nc=`head -n10 ${dfile} | gawk '/^ *[.0-9][.0-9]/ { printf "%d\n", NF-1; exit(1); }'`
echo "found ${nc} curves" 1>&2

# Obtain the max order:
npmax=`head -n 20 ${dfile} | egrep -e '^[#] *npmax *=' | sed -e 's:^.*[=] *::g'`
echo "npmax = ${npmax}" 1>&2 

# Create file with key frequencies:
ifile="${tmp}_i.txt"
rm -f ${ifile}
echo "${flo0} 2.0" >> ${ifile}
echo "${flo1} 2.0" >> ${ifile}
echo "${fhi1} 2.0" >> ${ifile}
echo "${fhi0} 2.0" >> ${ifile}
fnyq=`echo "${fsmp}/2.0" | bc -lq`  # Nyquist frequency.
echo "${fnyq} 2.0" >> ${ifile}

# Create the gnuplot command file:
gfile=${tmp}.gpl
sep="plot"
ic=0
rm -f ${gfile}
printf "%s"'\\\n' "${sep}" >> ${gfile}
printf "  \"${ifile}\" using 1:2 notitle with impulses lw 6 lc rgb '#ffccaa'" >> ${gfile}
sep=","
while [[ ${ic} -lt ${nc} ]]; do
  printf "%s"'\\\n' "${sep}" >> ${gfile}
  printf "  \"${dfile}\" using 1:%d" "$(( ${ic} + 2 ))" >> ${gfile}
  if [[ $(( ${ic} + 1 )) -eq ${nc} ]]; then
    ctit="np=${npmax}"
  else
    ctit="np=$(( ${ic} + 1 ))"
  fi
  printf "  title \"%s\"" "${ctit}" >> ${gfile}
  printf "  with lines lw 3" >> ${gfile}
  ic=$(( ${ic} + 1 ))
  sep=","
done
printf "\n" >> ${gfile}

tfile="${tmp}.png"
export GDFONTPATH=.:${HOME}/ttf
gnuplot <<EOF

set term png size 2800,1200 font "arial,24"
set output "${tfile}"
Gmin = 1.0e-9
Gmax = 1
set logscale y
fmin = 0.001
fmax = 1.000
set logscale x
set grid xtics
set grid ytics
set xrange [(0.9*fmin):(1.1*fmax)]
set yrange [(0.5*Gmin):(2.0*Gmax)]
load "${gfile}"
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
rm ${gfile}

