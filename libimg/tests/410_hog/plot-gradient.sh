#! /bin/bash
# Last edited on 2012-01-14 03:35:30 by stolfilocal

ffile="$1"; shift;

tfile=".tmp.png"
hfile="${ffile%%.*}-v.png"
echo "plotting ${ffile} (${nb} bins)"

gnuplot << EOF
  set term png large size 1800,1800
  set size ratio -1
  set output ".tmp.png"
  xgr(ix,iy) = ( abs(column(ix))+abs(column(iy)) == 0 ? 0/0 : 10*column(ix) )
  ygr(ix,iy) = ( abs(column(ix))+abs(column(iy)) == 0 ? 0/0 : 10*column(iy) )
  plot "${ffile}" using 1:2:(xgr(3,4)):(ygr(3,4)) title "hog" with vectors lw 2
  quit
EOF

echo "writing ${hfile}"
convert ${tfile} -resize '50%' ${hfile}
display ${hfile}

