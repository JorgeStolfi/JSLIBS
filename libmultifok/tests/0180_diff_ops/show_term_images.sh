#! /bin/bash
# Last edited on 2025-01-30 10:29:44 by stolfi

echo "=== show_terms_images.sh =============================" 1>&2
echo "$@" 1>&2

inFolder="$1"; shift;    # Directory where images are 
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
frameTag="`printf -- '-zf%08.4f-df%08.4f' ${zFoc} ${zDep}`"
outImagePrefix="`printf -- '%s-img-%s' ${outPrefix} ${runPrefix}`"

tmp="/tmp/$$"

# Images with photo, {Z} height average {zave}, {Z} height deviation {zdev},
# actual blurring indicator {shrp}, computed blurring {shrp}, and error {score-sharp}:

csimg=${inPrefix}${frameTag}-cs.ppm

avimg=${outImagePrefix}${frameTag}-av.pgm
gvimg=${outImagePrefix}${frameTag}-gv.pgm
dvimg=${outImagePrefix}${frameTag}-dv.pgm
nrimg=${outImagePrefix}${frameTag}-nr.pgm

# Omit ${csimg} for lack of space:
tmp1img="${tmp}-show1.ppm"
convert +append ${avimg} ${gvimg} ${dvimg} ${nrimg}  ${tmp1img}

azimg=${inPrefix}${frameTag}-az.pgm
dzimg=${inPrefix}${frameTag}-dz.pgm
shimg=${inPrefix}${frameTag}-sh.pgm
mkimg=${outImagePrefix}${frameTag}-mk.pgm

tmp2img="${tmp}-show2.pgm"
convert +append ${azimg} ${dzimg} ${shimg} ${mkimg} ${tmp2img}

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
