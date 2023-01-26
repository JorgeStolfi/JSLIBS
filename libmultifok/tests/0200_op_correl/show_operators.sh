#! /bin/bash
# Last edited on 2023-01-17 12:59:19 by stolfi

inPrefix="$1"; shift;
outPrefix="$1"; shift

tmp="/tmp/$$"

# Image with photo, relative {Z} height {zave}, {Z} height deviation {zdev},
# actual sharpness {sharp}, computed sharpness {score}, and error {score-sharp}:

cimg=${inPrefix}-c.ppm
simg=${inPrefix}-s.pgm
zimg=${inPrefix}-z.pgm
dimg=${inPrefix}-d.pgm
tczdimg="${tmp}-czd-w.ppm"
convert +append ${cimg} ${zimg} ${dimg} ${tczdimg}

fimg=${outPrefix}-f.pgm
eimg=${outPrefix}-e.pgm
tsfeimg="${tmp}-sfe-w.ppm"
convert +append ${simg} ${fimg} ${eimg} ${tsfeimg}

tczdsfeimg="${tmp}-czdsfex.ppm"
convert -append ${tczdimg} ${tsfeimg} ${tczdsfeimg}

# Images with photo, actual sharpness {sharp}, and basis coefficient {coeff[kb]}
for bimg in ${outPrefix}*-b.pgm ; do \
  tail="${bimg/${outPrefix}/}"
  tcsbimg="${tmp}${tail}-csby.ppm"
  echo "${bimg} --> ${tcsbimg}" 1>&2 
  convert +append ${cimg} ${simg} ${bimg} ${tcsbimg} 
done 

display -filter Box -resize '400%' ${tmp}*-*x.ppm ${tmp}*-*y.ppm
