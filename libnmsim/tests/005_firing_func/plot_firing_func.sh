#! /bin/bash
# Last edited on 2019-01-03 00:29:41 by jstolfi

fname="$1"; shift

tmp="/tmp/$$"

# END COMMAND LINE PARSING
# ----------------------------------------------------------------------

# Prefix for temporary file names
tmp="/tmp/$$"

dfile="${fname}.txt"
tfile="${tmp}.png"
pfile="${fname}.png"

export GDFONTPATH="."

gnuplot <<EOF
  set term png truecolor size 1600,800 font "arial,18"
  set output "${tfile}"
  set yrange [-0.001:1.001]
  set key left
  plot \
    "${dfile}" using 4:5 notitle with lines lt 1 lw 1 lc rgb '#ffcc77', \
    ""         using 1:2 title "Phi(V)" with lines lt 1 lw 1 lc rgb '#008800', \
    ""         using 1:(15*column(3)) title "15*Phi'(V)" with lines lt 1 lw 1 lc rgb '#cc0000'
  quit
EOF

if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile}
  display ${pfile}
  rm -f ${tfile}
fi
  
