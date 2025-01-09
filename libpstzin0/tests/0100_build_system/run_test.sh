#! /bin/bash
# Last edited on 2025-01-08 03:51:48 by stolfi

nmap="$1"; shift;
full="$1"; shift;

PROG="test_build_system"

size=0256x0256

indir="in/slope_to_height"

in_slopes_fni="${indir}/${nmap}-*-${size}-G.fni"
in_weights_fni="${indir}/${nmap}-*-${size}-W.fni"
outPrefix="out/${nmap}-${full}"
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
