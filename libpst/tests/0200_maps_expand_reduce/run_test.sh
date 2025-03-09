#! /bin/bash
# Last edited on 2025-02-28 05:24:59 by stolfi

fp="$1"; shift;

PROG="test_maps_shrink_expand"

make -k -f Makefile ${PROG}
  
printf "function ${fp} ..."
${PROG} -function ${fp/-*/}
for tag in GI GS-cmp GS-ref GS-dif ZI ZS-cmp ZS-ref ZS-dif ZE-cmp ZE-ref ZE-dif ; do
  ofiles=( out/${fp}-*-*-${tag}.fni )
  if [[ ${#ofiles[@]} -ne 1 ]]; then
    echo "bad ofiles = ( ${ofiles[*]} )" 1>&2; exit 1
  fi
  ofile="${ofiles[0]}"
  size="${ofile/-${tag}.fni}"
  size="${size/out\/${fp}-}"
  # echo "size = '${size}'" 1>&2
  nx="${size/-*}"; nx=$(( 10#${nx} + 0 ))
  ny="${size/*-}"; ny=$(( 10#${ny} + 0 ))
  if [[ ( "/${tag}" == "/GS-dif" ) || ( "/${tag}" == "/ZS-dif" ) ||  ( "/${tag}" == "/ZE-dif" ) ]]; then
    plot_cmp_ref.sh ${fp} ${size} ${tag/-dif*}
  else
    if [[ ( "/${tag}" == "/GI" ) || ( "/${tag}" == "/GS-ref" ) ||  ( "/${tag}" == "/GS-cmp" ) ]]; then
      scale=10
    else
      scale=`echo "sqrt(${nx}^2+${ny}^2)/2" | bc -lq`
      scale=`printf "%.2f" "${scale}"`
    fi
    fni_view -hist T -scale ${scale} ${ofile}
  fi
done
echo "OK"
