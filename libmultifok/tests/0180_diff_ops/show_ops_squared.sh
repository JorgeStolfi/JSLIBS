#! /bin/bash
# Last edited on 2023-01-25 23:43:31 by stolfi

inPrefix="$1"; shift;
outPrefix="$1"; shift

tmp="/tmp/$$"

# Image with photo, relative {Z} height {zave}, {Z} height deviation {zdev},
# actual sharpness {sharp}, computed sharpness {score}, and error {score-sharp}:

cimg=${inPrefix}-c.ppm
avimg=${outPrefix}-av.pgm
dvimg=${outPrefix}-dv.pgm
nimg=${outPrefix}-n.pgm
show1img="${tmp}-show1.ppm"
convert +append ${cimg} ${avimg} ${dvimg} ${nimg} ${show1img}

# Images with photo, actual sharpness {sharp}, and basis coefficient {coeff[kb]}
for bqimg in ${outPrefix}*-bq.pgm ; do \
  tail="${bqimg/${outPrefix}/}"
  show2img="${tmp}${tail}-show2.ppm"
  echo "${bqimg} --> ${show2img}" 1>&2 
  convert +append ${cimg} ${simg} ${bqimg} ${show2img} 
done 

display -filter Box -resize '400%' ${tmp}*-*show1.ppm ${tmp}*-*show2.ppm
