#! /usr/bin/bash
# Last edited on 2025-02-11 15:35:52 by stolfi

export GDFONTPATH="ttf"

gnuplot <<EOF
  set term postscript eps noenhanced color background '#ffffff' size 18,6 font "TimesRoman" 32
  set output "out/plot.eps"
  # set term X11 size 1800,600
  
  set xrange[-1.2:+1.2]
  set yrange[-2.1:+2.1]
  set samples 10000
  
  set zeroaxis
  
  max(x,y) = (x >= y ? x : y)
  clip(x) = (abs(x) > 1.0 ? 0.0 : x)
  goodnx(x) = -2*x
  goodnz(x) = sqrt(1-x**2)
  badnx(x) = -2*x/sqrt(1 + 4*x**2)
  badnz(x) = 1/sqrt(1 + 4*x**2)

  set multiplot layout 1,5
  set title "height"
  plot \
    1.0*sqrt(max(0,1-x**2)) notitle with lines lc rgb '#ff0000'

  set title "nx (true)"
  plot \
    1.0*goodnx(clip(x)) notitle with lines lc rgb '#ff0000'

  set title "nz (true)"
  plot \
    1.0*goodnz(clip(x)) notitle with lines lc rgb '#ff0000'

  set title "nx (wrong)"
  plot \
    1.0*badnx(clip(x)) notitle with lines lc rgb '#ff0000'

  set title "nz (wrong)"
  plot \
    1.0*badnz(clip(x)) notitle with lines lc rgb '#ff0000'

  unset multiplot
  pause mouse
EOF

if [[ -s out/plot.eps ]]; then
  evince out/plot.eps
fi
