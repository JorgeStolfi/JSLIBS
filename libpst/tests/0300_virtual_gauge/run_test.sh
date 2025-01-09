#! /bin/bash 
# Last edited on 2025-01-04 22:37:57 by stolfi

nf="$1"; shift

gauge_file="in/gauge-${nf}-parms.txt"
light_file="in/light-${nf}-parms.txt"

if [[ ! ( -s ${gauge_file} ) ]]; then echo "** file ${gauge_file} missing" 1>&2; exit 1; fi
if [[ ! ( -s ${light_file} ) ]]; then echo "** file ${light_file} missing" 1>&2; exit 1; fi

outPrefix="out/test-${nf}"
rm -f ${outPrefix}-*.{fni,png,ppm,pgm}

set -x
test_virtual_gauge \
  -gauge `cat ${gauge_file}` \
  -light `cat ${light_file}` \
  -outPrefix ${outPrefix}
set +x
  
sy_png_file="${outPrefix}-sy.png"
er_png_file="${outPrefix}-er.png"

if [[ -s ${sy_png_file} ]]; then
  display -title '%f' ${sy_png_file} 
else
  echo "** file ${sy_png_file} not created" 1>&2; exit 1
fi
if [[ -s ${er_png_file} ]]; then
  display -title '%f' ${er_png_file} 
else
  echo "** file ${er_png_file} not created" 1>&2; exit 1
fi
