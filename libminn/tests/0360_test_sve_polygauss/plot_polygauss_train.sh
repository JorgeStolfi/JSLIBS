#! /bin/bash
# Last edited on 2024-01-09 02:39:26 by stolfi

# Plots a polygaussian pulse and its component pulses.

show="$1"; shift
ifile="$1"; shift
np="$1"; shift

if [[ ! ( -s ${ifile} ) ]]; then echo "** file ${ifile} missing or empty" 1>&2; exit 1; fi
# Remove leading zeros (CROCK!)
np=`echo ${np} | sed -e 's:^0*::g'`
if [[ ${np} -lt 2 ]]; then echo "** invalid {np}" 1>&2; exit 1; fi


# Create the gnuplot command file:
color=( "ff0000" "996600" "338800" "008800" "007755" "0033ff" "5500ff" "aa0066" "aa2200" "557700" "117700" "007722" "005588" "2222ff" "880088" )
tmp="/tmp/$$"
gfile="${tmp}.gpl"
printf "plot %s\n" '\' > ${gfile}
ip=0
while [[ ${ip} -lt ${np} ]]; do
  printf "" >> ${gfile}
  printf "  \"%s\" using 2:%d title \"%02d\"" "${ifile}" $(( ${ip} + 4 )) ${ip} >> ${gfile}
  printf " with lines lw 2 lc rgb '#%s', %s\n" "${color[${ip}]}" '\' >> ${gfile}
  ip=$(( ${ip} + 1 ))
done
printf "  \"${ifile}\" using 2:3 title \"tot\" with lines lw 2 lc rgb '#000000'\n"  >> ${gfile}
  
for ysc in 0 1; do
  tfile="${tmp}.png"
  export GDFONTPATH=.:${HOME}/ttf
  gnuplot << EOF
  ysc = ${ysc} + 0
  set term png size 1600,1600 font "arial,24"
  set output "${tfile}"

  if (ysc != 0) {
    set yrange [-0.01:1.2]
    set xrange [-6.0:+7.0]
  } else {
    set yrange [0.97:1.03]
    set xrange [-1.0:+2.0]
  }

  unset logscale x
  unset logscale y
  
  set xtics
  set mxtics 5
  set grid xtics lt 1 lw 3 lc rgb '#dddddd', lt 1 lw 1.5 lc rgb '#dddddd'
  set grid mxtics
  
  set ytics
  set mytics 5
  set grid ytics lt 1 lw 3 lc rgb '#dddddd', lt 1 lw 1.5 lc rgb '#dddddd'
  set grid mytics
  
  load "${gfile}"
EOF
  if [[ -s ${tfile} ]]; then
    pfile="${ifile/.txt/}_${ysc}.png"
    convert ${tfile} -resize '50%' ${pfile}
    if [[ -s  ${pfile} ]]; then
      if [[ "/${show}" == "/YES" ]]; then
        display -title '%f' ${pfile}
      fi
      rm ${tfile}
    else
      echo "** convert ${tfile} to ${pfile} failed" 1>&2; exit 1
    fi
  else
    echo "** plot failed - ${tfile} empty or missing" 1>&2; exit 1
    rm -f ${tfile}
  fi
done
