#! /bin/bash
# Last edited on 2025-03-01 17:48:53 by stolfi

mapName="$1"; shift;
withHints="$1"; shift;

PROG="test_build_system"

indir="in/slope_to_height"

# 0001x0001 0032x0024 0256x0192
for size in 0032x0024 ; do

  nx="${size/x*/}"; nx=$(( 10#${nx} + 0 ))
  ny="${size/*x/}"; ny=$(( 10#${ny} + 0 ))
  rad=`echo "sqrt(${nx}*${nx}+${ny}*${ny})/2" | bc -lq`
  echo "nx = ${nx}  ny = ${ny}  rad = ${rad}" 1>&2

  in_slopes_fni="${indir}/${mapName}-${size}-G.fni"
  in_heights_fni="${indir}/${mapName}-${size}-Z.fni"

  outPrefix="out/${mapName}-${size}-hints${withHints}"

  out_system_txt="${outPrefix}-S.txt"
  out_weights_fni="${outPrefix}-SW.fni"

  if [[ "/${withHints}" == "/Y" ]]; then
    hintsOps=( "-hints" "${in_heights_fni}" )
  else
    hintsOps=()
  fi

  rm -f ${outPrefix}*.{txt,fni}
  set -x
  ${PROG} \
    -slopes ${in_slopes_fni} \
    ${hintsOps[@]} \
    -outPrefix ${outPrefix}
  set +x
  if [[ -s ${out_system_txt} ]]; then 
    ls -l ${out_system_txt} 1>&2
  else
    echo "** ${out_system_txt} not created" 1>&2 ; exit 1
  fi
  if [[ -s ${out_weights_fni} ]]; then 
    fni_view -hist T -scale ${rad} ${out_weights_fni}
  else
    echo "** ${out_weights_fni} not created" 1>&2 ; exit 1
  fi
done
echo "OK"
