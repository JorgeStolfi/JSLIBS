#! /bin/bash
# Last edited on 2025-01-14 15:44:44 by stolfi

mapName="$1"; shift;
full="$1"; shift;

PROG="test_build_system"

size=0256x0192

indir="in/slope_to_height"

in_slopes_fni="${indir}/${mapName}-${size}-G.fni"
in_weights_fni="${indir}/${mapName}-${size}-W.fni"
outPrefix="out/${mapName}-${size}-${full}"
out_system_txt="${outPrefix}-S.sys"
out_weights_pgm="${outPrefix}-W.pgm"

rm -f ${outPrefix}*.{sys,pgm}
set -x
${PROG} -slopes ${in_slopes_fni} -weights ${in_weights_fni} -outPrefix ${outPrefix} -full ${full}
set +x
if [[ -s ${out_system_txt} ]]; then 
  ls -l ${out_system_txt} 1>&2
else
  echo "** ${out_system_txt} not created" 1>&2 ; exit 1
fi
echo "OK"
