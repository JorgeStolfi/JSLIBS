#! /bin/bash
# Last edited on 2023-01-17 12:59:19 by stolfi

inPrefix="$1"; shift;
outPrefix="$1"; shift

tmp="/tmp/$$"

# Image with photo, relative {Z} height {zave}, {Z} height deviation {zdev},
# actual sharpness {sharp}, computed sharpness {score}, and error {score-sharp}:

csimg=${inPrefix}-cs.ppm
shimg=${inPrefix}-sh.pgm
azimg=${inPrefix}-az.pgm
dzimg=${inPrefix}-dz.pgm
tczdimg="${tmp}-czd-w.ppm"
convert +append ${csimg} ${azimg} ${dzimg} ${tczdimg}

scimg=${outPrefix}-sc.pgm
esimg=${outPrefix}-es.pgm
tsfeimg="${tmp}-sfe-w.ppm"
convert +append ${shimg} ${scimg} ${esimg} ${tsfeimg}

tczdsfeimg="${tmp}-czdsfex.ppm"
convert -append ${tczdimg} ${tsfeimg} ${tczdsfeimg}

# Images with photo, actual sharpness {sharp}, and basis coefficient {coeff[kb]}
for bcimg in ${outPrefix}*-bc.pgm ; do \
  tail="${bcimg/${outPrefix}/}"
  tcsbimg="${tmp}${tail}-csby.ppm"
  echo "${bcimg} --> ${tcsbimg}" 1>&2 
  convert +append ${csimg} ${shimg} ${bcimg} ${tcsbimg} 
done 

display -filter Box -resize '400%' ${tmp}*-*x.ppm ${tmp}*-*y.ppm
