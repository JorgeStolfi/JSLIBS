#! /bin/bash
# Last edited on 2023-03-07 14:06:50 by stolfi

ims=()
tts=()
for signed in 0 1 ; do
  if [[ ${signed} -gt 0 ]]; then sc='s'; ts='signed'; else sc='u'; ts='unsigned'; fi
  for style in 0 1 2 3 4 5 ; do
    for cycles in -5 -4 -3 -2 -1 0 +1 +2 +3 +4 +5 ; do 
      prefix="`printf 'out/path_%s%s_c%+03d' ${sc} ${style} ${cycles}`"
      txfile="${prefix}.txt"
      if [[ -s ${txfile} ]]; then 
        # There should be an image attached:
        imfile="${txfile/.txt/.ppm}"

        title="frgg_path_map_${ts}_${style}, cycles = ${cycles}"
        echo "${title}" 1>&2
        ../../plot_color_path.sh ${txfile} ${imfile} ${signed} "${title}"

        ims+=( ${imfile} )
        tts+=( "${ts}:${style}:${cycs}" )
      fi
    done
    ims+=( spacer.ppm )
    tts+=( "." )
  done
done

# Assemble all images in one

# convert -background "rgb(255,255,255)" +append ${ims[@]} out/.all.ppm
# display -rotate -90 out/.all.ppm
