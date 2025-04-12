#! /bin/bash
# Last edited on 2025-03-24 15:55:01 by stolfi

for ff in 0.9 1.a 2.b 3.c 5.e 8.k ; do
  f0="${ff/.*/}"; f1="${ff/*./}"
  echo "f0 = ${f0} f1 = ${f1}" 1>&2 
  convert projects/texture-bank/pgm-512x512/wavys-${f0}${f1}.pgm wavys${f0}${f1}.png 
    
  for m in 1 2 3 4 5 6 ; do
    echo "  m = ${m}" 1>&2 
    convert projects/texture-bank/pgm-512x512/spots-0${m}.pgm spots0${m}.png
    multiply_gray_textures.sh spots0${m} wavys${f0}${f1} pardo${m}${f0}
    display pardo${m}${f0}.png spots0${m}.png wavys${f0}${f1}.png
  done
    
done

