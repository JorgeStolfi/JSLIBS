#! /bin/bash
# Last edited on 2023-01-25 07:54:57 by stolfi

kt="$1"; shift

xkb=`printf "%03d" ${kt}`

export GDFONTPATH="${HOME}/tt-fonts"

gnuplot << EOF

set term X11 
kt = ${kt}

term(k) = column(k+5)
grad(dum) = term(2) + term(5)
lap(dum) = term(9) + term(14) + term(20)

set xrange [-0.01:+1.01]
plot "out/test-005.0000-bD--rdata.txt" using 4:(lap(0)) title "term[${kt}]" with points pt 7 lc rgb '#009900'


pause 300
EOF
