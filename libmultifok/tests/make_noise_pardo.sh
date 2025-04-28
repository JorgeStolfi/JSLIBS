#! /bin/bash
# Last edited on 2025-04-18 21:54:15 by stolfi

spot_sizes=( 1 2 3 4 5 6 )

for m in ${spot_sizes[@]} ; do
  echo "  m = ${m}" 1>&2 
  spots_bank="spots-0${m}"
  spots_here="spots0${m}"
  convert projects/texture-bank/pgm-512x512/${spots_bank}.pgm ${spots_here}.png

for ff in 0.9 1.a 2.b 3.c 5.e 8.k ; do
  f0="${ff/.*/}"; f1="${ff/*./}"
  echo "f0 = ${f0} f1 = ${f1}" 1>&2 
  noise_bank="wavys-${f0}${f1}"
  noise_here="wavys${f0}${f1}"
  convert projects/texture-bank/pgm-512x512/${noise_bank}.pgm ${noise_here}.png 
    
  for m in ${spot_sizes[@]} ; do
    echo "  m = ${m}" 1>&2 
    spots_here="spots0${m}"
    pardo_here="pardo${m}${f0}"
    multiply_gray_textures.sh ${spots_here} ${noise_here} ${pardo_here}
    display ${pardo_here}.png ${spots_here}.png ${noise_here}.png
  done
    
done

