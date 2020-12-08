#! /bin/bash
# Last edited on 2012-01-14 02:28:20 by stolfilocal

nb="$1"; shift;
dfile="$1"; shift;

tfile=".tmp.png"
hfile="${dfile%%.*}.png"
echo "plotting ${dfile} (${nb} bins)"

gnuplot << EOF
  set term png large size 1280,960
  set output ".tmp.png"
  set xrange [-0.1:2*pi+0.1]
  set yrange [-0.1:*]
  plot "${dfile}" using 2:5 title "hog" with histeps lw 2
  quit
EOF

echo "writing ${hfile}"
convert ${tfile} -resize '50%' ${hfile}
display ${hfile}

