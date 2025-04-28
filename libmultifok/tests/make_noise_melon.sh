#! /bin/bash
# Last edited on 2025-02-08 20:22:46 by stolfi

for ns in 2.512 3.256 4.128 5.64 6.32 ; do
  n="${ns/.*/}"; s="${ns/*./}"
  h=$(( ${s} / 2 ))
  convert noise01.png \
    -crop "${s}x${s}+${h}+${s}" \
    -resize '1024x1024' \
    -sigmoidal-contrast '80x50%' \
    noise0${n}.png
    
  display noise0${n}.png
  for m in 1 2 3 4 5 6 ; do
    if [[ ${m} -lt ${n} ]]; then
      multiply_gray_textures.sh noise0${m} noise0${n} melon${m}${n}
    fi
  done
    
done

