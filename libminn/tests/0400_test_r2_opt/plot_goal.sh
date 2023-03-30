#! /bin/bash
# Last edited on 2017-06-06 20:26:33 by stolfilocal

dfile="$1"; shift
dname="${dfile%.*}"

export GDFONTPATH=ttf

gnuplot <<EOF
  set term png medium size 1200,1200 font "courbd,18"
  set output "out/.tmp.png"
  set hidden3d
  splot "out/${dname}.dat" using 1:2:3 with lines 
EOF

if [[ -r out/.tmp.png ]]; then
  convert out/.tmp.png -resize '50%' out/${dname}.png 
  display -title '%f' out/${dname}.png 
else
  rm -f out/${dname}.png 
fi
