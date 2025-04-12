#! /bin/bash
# Last edited on 2025-03-19 15:17:34 by stolfi

# Plots a potential file created by {test_sve_charges}

nq="$1"; shift;
potfile="$1"; shift;

# Remove leading zeros from ${nq} (THE OCTAL CONVENTION IS A CROCK!!!!!!!)
nq=`echo "${nq}" | sed -e 's:^0*::'`

gnuplot << EOF
set terminal X11
nq=${nq}
set hidden3d
minpot = 1.07*(nq - 0.75)**2
# minpot = 200
set zrange [(0.80*minpot):(1.30*minpot)]
splot "${potfile}" using 1:2:3 title "potential" with lines lc rgb '#008800'
pause mouse
EOF
