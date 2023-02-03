#! /bin/bash
# Last edited on 2023-01-26 16:12:31 by stolfi

outPrefix="$1"; shift;    # Prefix for all images.
zDep_FMT="$1"; shift;  # Depth of focus, formatted as in image names.

tmp="/tmp/$$"

# Sharp images:

outPrefixSharp="${outPrefix}-sharp"
csimg="${outPrefixSharp}-cs.ppm"
azimg="${outPrefixSharp}-az.pgm"
dzimg="${outPrefixSharp}-dz.pgm"
showImgSharp="${outPrefixSharp}-show.ppm"
convert +append ${csimg} ${azimg} ${dzimg} ${showImgSharp}
if [[ ! ( -s ${showImgSharp} ) ]]; then
  echo "sharp montage failed" 1>&2; exit 1
fi

outPrefixStack="${outPrefix}-fd${zDep_FMT}"
for csimg in ${outPrefixStack}*-cs.ppm ; do
  outPrefixFrame="${csimg/-cs.ppm}"
  frameName="${outPrefixFrame/*zf/zf}"
  echo "=== ${frameName} ===" 1>&2

  shimg="${outPrefixFrame}-sh.pgm"
  azimg="${outPrefixFrame}-az.pgm"
  dzimg="${outPrefixFrame}-dz.pgm"
  showImgFrame="${outPrefixFrame}-show.ppm"
  
  tmpPrefixFrame="${tmp}-${frameName}"
  row1img="${tmpPrefixFrame}-cs.ppm"
  row2img="${tmpPrefixFrame}-zd.ppm"
  convert +append ${csimg} ${shimg} ${row1img}
  convert +append ${azimg} ${dzimg} ${row2img}
  convert -append ${row1img} ${row2img} ${showImgFrame}
  rm -f ${row1img} ${row2img}
  if [[ ! ( -s ${showImgFrame} ) ]]; then
    echo "sharp montage failed" 1>&2; exit 1
  fi
done
display -title '%f' -filter Box -resize '400%' ${showImgSharp} ${outPrefixStack}*-show.ppm
