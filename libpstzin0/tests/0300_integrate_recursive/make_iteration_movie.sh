#! /bin/bash
# Last edited on 2025-03-02 11:50:51 by stolfi

inPrefix="$1"; shift
reportStep="$1"; shift
tag="$1"; shift
vmin="$1"; shift
vmax="$1"; shift

tmp="/tmp/$$"

iterDir="${inPrefix}-iters"
files=( `( cd ${iterDir} && ls it-[0-9][0-9]-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]-${tag}.fni ) | sort -t- -k2nr -k3n` )
frames=()
for f in ${files[@]}; do
  # echo ${f} 1>&2
  levit="${f/-${tag}.fni/}"; levit="${levit/it-/}";
  level="${levit/-*/}"; iter="${levit/*-/}";
  echo "level = '${level}'  iter = '${iter}'" 1>&2
  niter=$(( 10#${iter} ))
  rem=$(( ${niter} % ${reportStep} ))
  echo "niter = '${niter}'  rem = '${rem}'" 1>&2
  if [[ ${rem} -eq 0 ]]; then
    echo "plotting $f ..." 1>&2
    pfile="${tmp}-${num}.png"
    fni_plot.sh -channel 0 -range ${vmin} ${vmax} -title "iteration ${num}" < ${iterDir}/$f > ${pfile}
    frames+=( ${pfile} )
  fi
done
gfile="${tag}.gif"
rm -fv ${iterDir}/${gfile}
echo "creating ${gfile} ..." 1>&2
( cd ${iterDir} && convert -delay 50 ${frames[0]} ${frames[0]} ${frames[0]} ${frames[@]} ${frames[-1]} ${frames[-1]} ${frames[-1]} ${gfile} )
rm ${tmp}-*.png
