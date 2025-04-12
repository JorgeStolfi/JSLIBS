#! /bin/bash
# Last edited on 2025-03-22 06:08:12 by stolfi

outDir="$1"; shift
imgName="$1"; shift

frames=( ` ( cd ${outDir} && ls ${imgName}-[0-9][0-9][0-9][0-9][0-9].png ) | sort ; echo ${imgName}-final.png ` ) 

outFile="${outDir}/${imgName}.mp4"

ffmpeg \
  -framerate 10 \
  -start_number 0 \
  -i "${outDir}/${imgName}-%05d.png" \
  -y \
  -r 10 \
  -vcodec libx264 \
  -crf 25 \
  -pix_fmt yuv420p \
  ${outFile}

mplayer ${outFile}
