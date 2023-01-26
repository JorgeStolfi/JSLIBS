#! /bin/bash
# Last edited on 2023-01-25 19:26:11 by stolfi

prefix="$1"; shift;

tmp="/tmp/$$"
mkdir ${tmp}

for cimg in ${prefix}*-c.ppm ; do
  simg=${cimg/-c.ppm/-s.pgm}
  zimg=${cimg/-c.ppm/-z.pgm}
  dimg=${cimg/-c.ppm/-d.pgm}
  
  tpref="${cimg/out\/*\//${tmp}}"
  tpref="${tpref/-c.ppm/}"
  
  echo "=== ${tpref} ===" 1>&2
  
  tcsimg=${tpref}-cs.ppm
  tzdimg=${tpref}-zd.ppm
  showimg=${tpref}-cszd.ppm
  convert +append ${cimg} ${simg} ${tcsimg}
  convert +append ${zimg} ${dimg} ${tzdimg}
  convert -append ${tcsimg} ${tzdimg} ${showimg}
  rm -f ${tcsimg} ${tzdimg}
done
display -title '%f' -filter Box -resize '400%' ${tpref}*-cszd.ppm
