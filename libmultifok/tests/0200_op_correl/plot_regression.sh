#! /bin/bash
# Last edited on 2023-01-09 03:42:56 by stolfi

fname="$1"; shift  # File name minus the "-b{B}-vals.txt" tail.
btype="$1"; shift  # "N", "D", or "H".
nw="$1"; shift     # Window size.

if [[ ${nw} -ne 3 ]]; then echo "** window size must be 3" 1>&2; exit 1; fi
nb=$(( $nw * $nw ))

tmp="/tmp/$$"

rfile="${fname}-b${btype}-regr.txt"  # Regression result file.
ffile="${fname}-b${btype}-form.txt"  # Formula file.

export GDFONTPATH="${HOME}/tt-fonts"

gnuplot << EOF

set term X11 

set xrange [-0.01:+1.01]
plot "${rfile}" using 2:3 title "fitted" with points pt 7 lc rgb '#009900'


pause 300
EOF
