#! /bin/bash
# Last edited on 2025-03-07 07:56:52 by stolfi

scene="$1"; shift;
size="$1"; shift;

dtseq_G="$1"; shift;
method_G="$1"; shift;

dtseq_H="$1"; shift;
method_H="$1"; shift;

dtseq_R="$1"; shift;
method_R="$1"; shift;

initial_opt="$1"; shift;  # "zero", "hints", or "reference"
initial_noise="$1"; shift; 
maxLevel="$1"; shift

PROG="test_integrate_recursive"

( cd ~/programs/c/JSLIBS/libpst && make -k -f Makefile uninstall build-lib install )
( cd ~/programs/c/JSLIBS/libpstzin0 && make -k -f Makefile uninstall build-lib install )
rm -f ${PROG}
make -k -f Makefile ${PROG}
if [[ ! ( -s ${PROG} ) ]]; then "** compilation errors - aborted" 1>&2; exit 1 ; fi

inDir="in/gustavo"

tmp="/tmp/$$"

testName="${dtseq_G}-${scene}-${method_G}-${size}"

in_slopes_fni="${inDir}/${testName}-G.fni"
in_normals_fni="${inDir}/${testName}-N.fni"
if [[ -s ${in_slopes_fni} ]]; then 
  echo "using slope map {G} = ${in_slopes_fni}}"  1>&2
  slonorm_opts=( -slopes ${in_slopes_fni} scale +10 +10 )
  input_tgfi='input G':"${in_slopes_fni}"
elif [[ -s ${in_normals_fni} ]]; then 
  echo "using normal map {N} = ${in_normals_fni}}"  1>&2
  slonorm_opts=( -normals ${in_normals_fni} scale +10 -10 )
  input_tgfi='input N':"${in_normals_fni}"
else
  echo "** can find neither ${in_slopes_fni} nor ${in_normals_fni}" 1>&2 ; exit 1
fi

in_hints_fni="${inDir}/${dtseq_H}-${scene}-${method_H}-${size}-Z.fni"
if [[ -s ${in_hints_fni} ]]; then
  echo "using hints map {H} = ${in_hints_fni}}"  1>&2
  hints_opts=( -hints ${in_hints_fni} scale 192 0.0 )
else
  echo "!! cannot find hints file ${in_hints_fni}" 1>&2
  hints_opts=( )
fi

in_refz_fni="${inDir}/${dtseq_R}-${scene}-${method_R}-${size}-Z.fni"
if [[ -s ${in_refz_fni} ]]; then
  echo "using reference map {R} = ${in_refz_fni}}"  1>&2
  reference_opts=( -reference ${in_refz_fni} scale 192 )
else
  reference_opts=( )
  echo "!! cannot find reference file ${in_refz_fni}" 1>&2
fi 

outPrefix="out/${testName}"
iterDir="${outPrefix}-iters"
iterPrefix="${iterDir}/it"

if [[ ${size} == "0512x0422" ]]; then
  clear_opts=( -clear 0 511 0 40 )
else
  clear_opts=( )
fi

rm -f ${outPrefix}*.{fni,sys,pgm,txt,png} 
rm -rf ${iterDir}
set -x
${PROG} \
  ${slonorm_opts[@]} \
  ${hints_opts[@]}  \
  ${reference_opts[@]} \
  ${clear_opts[@]} \
  -maxLevel ${maxLevel} \
  -initial ${initial_opt} ${initial_noise} \
  -outPrefix ${outPrefix} \
  -maxIter 500 \
  -reportStep 10 \
  -verbose
status=$?
set +x

if [[ ${status} -ne 0 ]]; then
  echo "** command failed with status=${status}" 1>&2; exit 1
fi

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
    echo "** ${fni_file} not found" 1>&2
  fi
}

echo "displaying the inputs ..." 1>&2
for tgfi in "${input_tgfi}" 'input H':"${in_hints_fni}" 'input R':"${in_refz_fni}" ; do
  showfni "${tgfi}"
done

echo "displaying the intermediate states ..." 1>&2
for level in 20 19 18 17 16 15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00 ; do
  levelPrefix="${outPrefix}-${level}"
  for stg in 'beg-G' 'beg-H' 'beg-Z' 'beg-E' 'end-Z' 'end-E' ; do
    fni_file="${levelPrefix}-${stg}.fni"
    xstg="${stg/-/ }"
    showfni "${xstg}:${fni_file}"
  done
done

exit 0

make_iteration_movie.sh ${outPrefix} 10 Z -${zscale} +${zscale}
make_iteration_movie.sh ${outPrefix} 10 E -2.0 +2.0
echo "OK"

