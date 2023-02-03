#! /bin/bash
# Last edited on 2023-01-30 18:22:59 by stolfi

echo "=== show_terms_images.sh =============================" 1>&2
echo "$@" 1>&2

inDir="$1"; shift;    # Directory where images are 
outPrefix="$1"; shift

sceneType="$1"; shift;
NX="$1"; shift;
NY="$1"; shift;
pattern="$1"; shift;
zDep="$1"; shift;
zFoc="$1"; shift;

sizeTag="`printf '%04dx%04d' ${NX} ${NY}`"
runPrefix="st${sceneType}-${sizeTag}-${pattern}"
inPrefix="in/${runPrefix}/frame"
frameTag="`printf -- '-fd%05.2f-zf%05.2f' ${zDep} ${zFoc}`"
outImagePrefix="`printf -- '%s-img-%s' ${outPrefix} ${runPrefix}`"

tmp="/tmp/$$"

# Image with photo, relative {Z} height {zave}, {Z} height deviation {zdev},
# actual sharpness {sharp}, computed sharpness {score}, and error {score-sharp}:

csimg=${inPrefix}${frameTag}-cs.ppm
azimg=${inPrefix}${frameTag}-az.pgm
dzimg=${inPrefix}${frameTag}-dz.pgm
shimg=${inPrefix}${frameTag}-sh.pgm
avimg=${outImagePrefix}${frameTag}-av.pgm
dvimg=${outImagePrefix}${frameTag}-dv.pgm
nrimg=${outImagePrefix}${frameTag}-nr.pgm
mkimg=${outImagePrefix}${frameTag}-mk.pgm

tmp1img="${tmp}-show1.ppm"
convert +append ${csimg} ${shimg} ${azimg} ${dzimg} ${tmp1img}

tmp2img="${tmp}-show2.pgm"
convert +append ${avimg} ${dvimg} ${nrimg} ${mkimg} ${tmp2img}

showimg="${outImagePrefix}${frameTag}-show.ppm"
convert -append ${tmp1img} ${tmp2img} ${showimg}

# Images with photo, actual sharpness {sharp}, and quadtratic term {term[kt]}
for tmimg in ${outImagePrefix}*-tm.pgm ; do \
  tail="${tmimg/*-kt/-kt}"
  tail="${tail/.pgm/}"
  tmp3img="${tmp}${tail}-show3.ppm"
  echo "${tmimg} --> ${tmp3img}" 1>&2 
  convert +append ${csimg} ${shimg} ${tmimg} ${tmp3img} 
done 

display -title '%f' -filter Box -resize '400%' ${showimg} `ls ${tmp}*-*show3.ppm | sort`
