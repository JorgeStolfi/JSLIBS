#! /bin/bash
# Last edited on 2025-03-03 03:10:43 by stolfi

func="$1"; shift;
size="$1"; shift;
hints="$1"; shift;
reference="$1"; shift;  # 0 
initial_opt="$1"; shift;  # "zero", "hints", or "reference"
initial_noise="$1"; shift; 
maxLevel="$1"; shift

PROG="test_integrate_recursive"

inDir="in/slope_to_height"

testName="${func}-${size}"

tmp="/tmp/$$"

in_slopes_fni="${inDir}/${testName}-G.fni"
if [[ ${reference} -ne 0 ]]; then
  in_refz_fni="${inDir}/${testName}-Z.fni"
  out_errz_fni="${outPrefix}-00-end-E.fni"
  reference_ops=( "-reference" ${in_refz_fni} )
else
  in_refz_fni="NONE/NONE"
  out_errz_fni="NONE/NONE"
  reference_ops=( )
fi
if [[ ${hints} -ne 0 ]]; then
  in_hints_fni="${inDir}/${testName}-Z.fni"
  hints_ops=( "-hints" ${in_hints_fni} 0.01 )
else
  hints_ops=()
fi

outPrefix="out/${testName}"
iterPrefix="${outPrefix}-iters/it"

rm -f {${outPrefix},${iterPrefix}}*.{fni,sys,pgm,txt,png}
set -x
${PROG} \
  -slopes ${in_slopes_fni} \
  ${hints_ops[@]}  \
  ${reference_ops[@]} \
  -maxLevel ${maxLevel} \
  -initial ${initial_opt} ${initial_noise} \
  -outPrefix ${outPrefix} \
  -reportStep 10 \
  -verbose
set +x

function showfni() {
  tgfi="$1"; shift

  tag="${tgfi/:*}"
  fni_file="${tgfi/*:}"
  fni_name="${fni_file/*\//}"
  
  if [[ -s ${fni_file} ]]; then
    echo "showing ${fni_file} ..." 1>&2
    ny="`cat ${fni_file} | egrep -e 'NY *='`"
    ny="${ny/*=/}"; ny="${ny/* /}"; 
    ny=$(( 10#${ny} + 0 ))

    if [[ "/${tag/* /}" == "/G"  ]]; then
      scale=1;
    else 
      scale=$(( ${ny} / 2 ));
    fi
    echo "tag = '${tag}'  ny = ${ny}  scale=${scale}" 1>&2
    pgm_file="${tmp}.pgm"
    png_file="${fni_file/.fni/.png}"
    fni_view -title "${tag} : ${fni_name}" -hist y -scale ${scale} ${fni_file}
    # fni_to_pnm -channel 0 -yAxis up < ${fni_file} > ${pgm_file}
    # convert ${pgm_file} ${png_file}
    # display -title "${tag} : ${fni_name}" -filter box -resize 'x800<' ${png_file}
  else
    echo "** ${fni_file} not created" 1>&2
  fi
}

for tgfi in 'input G':${in_slopes_fni} 'input H':${in_hints_fni} 'input R':${in_refz_fni} ; do
  showfni "${tgfi}"
done

for level in 20 19 18 17 16 15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00 ; do
  levelPrefix="${outPrefix}-${level}"
  for state in "beg" "end"; do 
    for tag in 'Z' 'G' 'E' ; do
      if [[ ( "/${state}" == "/beg" ) || ( "/${tag}" != "/G" ) ]]; then
        fni_file="${levelPrefix}-${state}-${tag}.fni"
        stg="${state} ${tag}"
        showfni "${stg}:${fni_file}"
      fi
    done
  done
done

exit 0

make_iteration_movie.sh ${outPrefix} 10 Z -${zscale} +${zscale}
make_iteration_movie.sh ${outPrefix} 10 E -2.0 +2.0
echo "OK"

