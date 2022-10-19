#! /bin/bash
# Last edited on 2021-12-19 07:32:12 by stolfi

datafile="$1"; shift

plotfile="${datafile/.dat/.png}"

tdatafile="out/.tmp.dat"
tplotfile="out/.tmp.png"

rm -f ${tdatafile}

cat ${datafile} | sort -k1,1g -k2,2g | break_3D_plot_lines.gawk -v field=1 > ${tdatafile}
echo "" >> ${tdatafile}
cat ${datafile} | sort -k2,2g -k1,1g | break_3D_plot_lines.gawk -v field=2 >> ${tdatafile}

gnuplot <<EOF
  set term pngcairo size 1200,1200
  set output "${tplotfile}"
  set xrange [-5:+5]
  set yrange [-5:+5]
  set zrange [-0.001:]
  splot "${tdatafile}" using 1:2:3 with lines, \
        "${tdatafile}" using 1:2:(0) with lines
EOF

if [[ -r ${tplotfile} ]]; then
  convert ${tplotfile} -resize '50%' ${plotfile} 
  display -title '%f' ${plotfile}
fi
