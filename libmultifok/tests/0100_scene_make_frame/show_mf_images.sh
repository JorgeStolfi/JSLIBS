#! /bin/bash
# Last edited on 2025-02-01 16:24:29 by stolfi

stackFolder="$1"; shift;    # Prefix for all images.
zDep_FMT="$1"; shift;    # Depth of focus, formatted as in image names.

tmp="/tmp/$$"

# ----------------------------------------------------------------------
# Assemble the sharp images:

tmpPrefixSharp="${tmp}-sharp"

sVal_sharp="${stackFolder}/frame-sharp/sVal.png"
iDeb="${stackFolder}/selected-pixels.png"

row1img="${tmpPrefixSharp}-sVal-iDeb.png"
convert +append ${sVal_sharp} ${iDeb} ${row1img}

hAvg_sharp="${stackFolder}/frame-sharp/hAvg.png"
hDev_sharp="${stackFolder}/frame-sharp/hDev.png"
sNrm_sharp="${stackFolder}/frame-sharp/sNrm.png"

row2img="${tmpPrefixSharp}-hAvg-hDev-sNrm.png"
convert +append ${hAvg_sharp} ${hDev_sharp} ${sNrm_sharp} ${row2img}

showImgSharp="${stackFolder}/frame-sharp/show.png"
convert -append ${row1img} ${row2img} ${showImgSharp}
if [[ ! ( -s ${showImgSharp} ) ]]; then
  echo "sharp montage failed" 1>&2; exit 1
fi

# ----------------------------------------------------------------------
# Assemble the per-frame images:

frameDirs=( `cd ${stackFolder} && ls -d frame-*-df${zDep_FMT} | sort` )
for frameFolder in ${frameDirs[@]} ; do
  tmpPrefixFrame="${tmp}-${frameFolder}"
  echo "=== ${frameFolder} ===" 1>&2

  sVal="${stackFolder}/${frameFolder}/sVal.png"
  shrp="${stackFolder}/${frameFolder}/shrp.png"

  row1img="${tmpPrefixFrame}-sVal-shrp.png"
  convert +append ${sVal} ${shrp} ${row1img}

  hAvg="${stackFolder}/${frameFolder}/hAvg.png"
  hDev="${stackFolder}/${frameFolder}/hDev.png"

  sNrm="${stackFolder}/${frameFolder}/sNrm.png"

  row2img="${tmpPrefixFrame}-hAvg-hDev-sNrm.png"
  convert +append ${hAvg} ${hDev} ${sNrm} ${row2img}

  showImgFrame="${stackFolder}/${frameFolder}/show.png"
  convert -append ${row1img} ${row2img} ${showImgFrame}
  rm -f ${row1img} ${row2img}
  if [[ ! ( -s ${showImgFrame} ) ]]; then
    echo "frame ${frameFolder} montage failed" 1>&2; exit 1
  fi
done

# ----------------------------------------------------------------------
# Show the assembled images:

( cd ${stackFolder} && display -geometry '1600x800' -title '%d/%f' -filter Box -resize '400%' frame-{sharp,*-df${zDep_FMT}}/show.png )
