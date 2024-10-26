#! /bin/bash
# Last edited on 2024-10-22 10:57:40 by stolfi

stackDir="$1"; shift;    # Prefix for all images.
zDep_FMT="$1"; shift;    # Depth of focus, formatted as in image names.

tmp="/tmp/$$"

# ----------------------------------------------------------------------
# Assemble the sharp images:

tmpPrefixSharp="${tmp}-sharp"

sVal_sharp="${stackDir}/frame-sharp/sVal.png"
iDeb="${stackDir}/selected-pixels.png"

row1img="${tmpPrefixSharp}-sVal-iDeb.png"
convert +append ${sVal_sharp} ${iDeb} ${row1img}

hAvg_sharp="${stackDir}/frame-sharp/hAvg.png"
hDev_sharp="${stackDir}/frame-sharp/hDev.png"

row2img="${tmpPrefixSharp}-hAvg-hDev.png"
convert +append ${hAvg_sharp} ${hDev_sharp} ${row2img}

showImgSharp="${stackDir}/frame-sharp/show.png"
convert -append ${row1img} ${row2img} ${showImgSharp}
if [[ ! ( -s ${showImgSharp} ) ]]; then
  echo "sharp montage failed" 1>&2; exit 1
fi

# ----------------------------------------------------------------------
# Assemble the per-frame images:

frameDirs=( `cd ${stackDir} && ls -d frame-*-df${zDep_FMT} | sort` )
for frameDir in ${frameDirs[@]} ; do
  tmpPrefixFrame="${tmp}-${frameDir}"
  echo "=== ${frameDir} ===" 1>&2

  sVal="${stackDir}/${frameDir}/sVal.png"
  shrp="${stackDir}/${frameDir}/shrp.png"

  row1img="${tmpPrefixFrame}-sVal-shrp.png"
  convert +append ${sVal} ${shrp} ${row1img}

  hAvg="${stackDir}/${frameDir}/hAvg.png"
  hDev="${stackDir}/${frameDir}/hDev.png"

  row2img="${tmpPrefixFrame}-hAvg-hDev.png"
  convert +append ${hAvg} ${hDev} ${row2img}

  showImgFrame="${stackDir}/${frameDir}/show.png"
  convert -append ${row1img} ${row2img} ${showImgFrame}
  rm -f ${row1img} ${row2img}
  if [[ ! ( -s ${showImgFrame} ) ]]; then
    echo "frame ${frameDir} montage failed" 1>&2; exit 1
  fi
done

# ----------------------------------------------------------------------
# Show the assembled images:

( cd ${stackDir} && display -geometry '1600x800' -title '%d/%f' -filter Box -resize '400%' frame-{sharp,*-df${zDep_FMT}}/show.png )
