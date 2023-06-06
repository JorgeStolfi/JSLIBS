#! /bin/bash
# Last edited on 2023-04-28 17:47:52 by stolfi

outPrefix="$1"; shift;    # Prefix for all images.
zDep_FMT="$1"; shift;  # Depth of focus, formatted as in image names.

tmp="/tmp/$$"

# ----------------------------------------------------------------------
# Assemble the sharp images:

outPrefixSharp="${outPrefix}-sharp"
tmpPrefixSharp="${tmp}-sharp"

csimg_sharp="${outPrefixSharp}-cs.ppm"
psimg_sharp="${outPrefixSharp}-ps.ppm"

row1img="${tmpPrefixSharp}-cs-ps.ppm"
convert +append ${csimg_sharp} ${psimg_sharp} ${row1img}

azimg_sharp="${outPrefixSharp}-az.pgm"
dzimg_sharp="${outPrefixSharp}-dz.pgm"

row2img="${tmpPrefixSharp}-az-dz.ppm"
convert +append ${azimg_sharp} ${dzimg_sharp} ${row2img}

showImgSharp="${outPrefixSharp}-show.ppm"
convert -append ${row1img} ${row2img} ${showImgSharp}
if [[ ! ( -s ${showImgSharp} ) ]]; then
  echo "sharp montage failed" 1>&2; exit 1
fi

# ----------------------------------------------------------------------
# Assemble the per-frame images:

outPrefixStack="${outPrefix}-fd${zDep_FMT}"
for csimg in ${outPrefixStack}*-cs.ppm ; do
  outPrefixFrame="${csimg/-cs.ppm}"
  frameName="${outPrefixFrame/*zf/zf}"
  tmpPrefixFrame="${tmp}-${frameName}"
  echo "=== ${frameName} ===" 1>&2

  shimg="${outPrefixFrame}-sh.pgm"

  row1img="${tmpPrefixFrame}-cs-sh.ppm"
  convert +append ${csimg} ${shimg} ${row1img}

  azimg="${outPrefixFrame}-az.pgm"
  dzimg="${outPrefixFrame}-dz.pgm"

  row2img="${tmpPrefixFrame}-az-dz.ppm"
  convert +append ${azimg} ${dzimg} ${row2img}

  showImgFrame="${outPrefixFrame}-show.ppm"
  convert -append ${row1img} ${row2img} ${showImgFrame}
  rm -f ${row1img} ${row2img}
  if [[ ! ( -s ${showImgFrame} ) ]]; then
    echo "sharp montage failed" 1>&2; exit 1
  fi
done

# ----------------------------------------------------------------------
# Show the assembled images:

display -title '%f' -filter Box -resize '400%' ${showImgSharp} ${outPrefixStack}*-show.ppm
