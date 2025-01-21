#! /bin/bash
# Last edited on 2025-01-16 08:35:26 by stolfi

mapName="$1"; shift;
keepNull="$1"; shift;
size="$1"; shift;

PROG="test_integrate_iterative"

inDir="in/slope_to_height"

testName="${mapName}-${size}"

in_slopes_fni="${inDir}/${testName}-G.fni"
in_weights_fni="${inDir}/${testName}-W.fni"
in_refz_fni="${inDir}/${testName}-Z.fni"

outPrefix="out/${testName}-${keepNull}"
out_heights_fni="${outPrefix}-dbg-Z.fni"
out_errz_fni="${outPrefix}-dbg-eZ.fni"
out_U_fni="${outPrefix}-dbg-U.fni"

rm -f ${outPrefix}*.{fni,sys,pgm,txt,png}
set -x
${PROG} \
  -slopes ${in_slopes_fni} \
  -weights ${in_weights_fni} \
  -compareZ ${in_refz_fni} \
  -outPrefix ${outPrefix} \
  -keepNull ${keepNull} \
  -reportStep 10 \
  -verbose
set +x
for ofile in ${out_heights_fni} ${out_heights_fni} ${out_U_fni} ; do
  if [[ -s ${ofile} ]]; then 
    if [[ "/${ofile}" == "/${out_errz_fni}" ]]; then scale=1000; else scale=1; fi
    fni_view -scale ${scale} ${ofile}
    pfile="${ofile/.fni/.pgm}"
    fni_to_pnm -yAxis up < ${ofile} > ${pfile}
    display -title '%f' -filter box -resize 'x800<' ${pfile}
  else
    echo "** ${ofile} not created" 1>&2
  fi
done
echo "OK"

