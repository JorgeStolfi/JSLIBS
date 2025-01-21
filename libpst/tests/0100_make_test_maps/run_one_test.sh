#! /bin/bash
# Last edited on 2025-01-19 23:43:17 by stolfi

func_num="$1"; shift
func_name="$1"; shift
nx="$1"; shift
ny="$1"; shift
noisy="$1"; shift  # "Y" or "N".

echo "creating images of function ${func_num} = '${func_name}' size ${nx} x ${ny}"

PROGDIR="."
PROG="make_test_maps"

if [[ "/${noisy}" == "/Y" ]]; then 
  noiseG=0.300 
elif [[ "/${noisy}" == "/N" ]]; then 
  noiseG=0.000
else
  echo "**invalid noisy = '${noisy}'" 1>&2; exit 1
fi
noiseW=0.000

if [[ ( ${nx} -le 1 ) || ( ${ny} -le 1 ) ]]; then smoothG=( 0 1 ); else smoothG=( 1  5 ); fi
if [[ ( ${nx} -le 1 ) || ( ${ny} -le 1 ) ]]; then smoothZ=( 1 5 ); else smoothZ=( 2  11 ); fi

outPrefix="out/${func_num}-${func_name}-${noisy}-"`printf "%04dx%04d" ${nx} ${ny}`
rm -f ${outPrefix}-*.{fni,ppm,pgm,eps}
set -x
${PROGDIR}/${PROG} \
  -function ${func_num} \
  -size ${nx} ${ny} \
  -numGrad Y \
  -maxGrad 100.0 \
  -maxGDiff 1.000 \
  -noiseG ${noiseG} \
  -outPrefix ${outPrefix}
set +x

ls -d ${outPrefix}-*.fni

out_heights_fni="${outPrefix}-Z.fni"
out_slopes_fni="${outPrefix}-G.fni"
out_normals_fni="${outPrefix}-N.fni"

view_fni=1
make_pnm=1
make_hist=0
make_plot=0

pfiles=()
efiles=()
hfiles=()

scale=`echo "${ny}/2" | bc -lq`

if [[ ${view_fni} -gt 0 ]]; then
  ofile=${out_heights_fni}; if [[ -s ${ofile} ]]; then fni_view -channel 0 -scale ${scale} ${ofile}; fi
  ofile=${out_slopes_fni};  if [[ -s ${ofile} ]]; then fni_view -channel 0 -scale ${scale} ${ofile}; fi
  ofile=${out_normals_fni}; if [[ -s ${ofile} ]]; then fni_view -channel 0 -scale ${scale} -texture ${ofile} ${ofile}; fi
fi

if [[ ${make_pnm} -gt 0 ]]; then
  ofile=${out_heights_fni}; if [[ -s ${ofile} ]]; then 
    pfile="${ofile/.fni/.pgm}"
    fni_to_pnm -channel 0 < ${ofile} > ${pfile}
    if [[ -s ${pfile} ]]; then pfiles+=( ${pfile} ); fi
    wfile="${ofile/.fni/W.pgm}"
    fni_to_pnm -channel 1 < ${ofile} > ${wfile}
    if [[ -s ${wfile} ]]; then pfiles+=( ${wfile} ); fi
  fi
  
  ofile=${out_slopes_fni};  if [[ -s ${ofile} ]]; then 
    pfile="${ofile/.fni/.ppm}"
    fni_to_pnm -channels 0 1 0 -center 0 -uniform < ${ofile} > ${pfile}
    if [[ -s ${pfile} ]]; then pfiles+=( ${pfile} ); fi
    wfile="${ofile/.fni/W.pgm}"
    fni_to_pnm -channel 2 < ${ofile} > ${wfile}
    if [[ -s ${wfile} ]]; then pfiles+=( ${wfile} ); fi
  fi
  
  ofile=${out_normals_fni}; if [[ -s ${ofile} ]]; then 
    pfile="${ofile/.fni/.ppm}"
    fni_to_pnm -channels 0 1 2 -min -1 -max +1 < ${ofile} > ${pfile}
    if [[ -s ${pfile} ]]; then pfiles+=( ${pfile} ); fi
    wfile="${ofile/.fni/W.pgm}"
    fni_to_pnm -channel 3 < ${ofile} > ${wfile}
    if [[ -s ${wfile} ]]; then pfiles+=( ${wfile} ); fi
  fi
fi  

if [[ ${make_plot} -gt 0 ]]; then

  ofile=${out_heights_fni}; if [[ -s ${ofile} ]]; then 
    efile="${ofile/.fni/.eps}"; efiles+=( ${efile} )
    fni_plot.sh -channel 0 -title "Z" < ${ofile} > ${efile}
    if [[ -s ${efile} ]] ; then efiles+=( ${efile} ); fi
  fi

  ofile=${out_slopes_fni};  if [[ -s ${ofile} ]]; then 
    efile_X="${ofile/.fni/X.eps}"; efiles+=( ${efile_X} )
    fni_plot.sh -channel 0 -title "dZ/dX" < ${ofile} > ${efile_X}
    if [[ -s ${efile_X} ]] ; then efiles+=( ${efile_X} ); fi

    efile_Y="${ofile/.fni/Y.eps}"; efiles+=( ${efile_Y} )
    fni_plot.sh -channel 1 -title "dZ/dY" < ${ofile} > ${efile_Y}
    if [[ -s ${efile_Y} ]] ; then efiles+=( ${efile_Y} ); fi
  fi

  ofile=${out_normals_fni}; if [[ -s ${ofile} ]]; then 
    efile_X="${ofile/.fni/X.eps}"; 
    fni_plot.sh -channel 0 -title "nrm.X" < ${ofile} > ${efile_X};
    if [[ -s ${efile_X} ]] ; then efiles+=( ${efile_X} ); fi

    efile_Y="${ofile/.fni/Y.eps}";
    fni_plot.sh -channel 1 -title "nrm.Y" < ${ofile} > ${efile_Y}
    if [[ -s ${efile_Y} ]] ; then efiles+=( ${efile_Y} ); fi

    efile_Z="${ofile/.fni/Z.eps}"; 
    fni_plot.sh -channel 2 -title "nrm.Z" < ${ofile} > ${efile_Z}
    if [[ -s ${efile_Z} ]] ; then efiles+=( ${efile_Z} ); fi
  fi

fi

if [[ ${make_hist} -gt 0 ]]; then

  ofile=${out_heights_fni}; if [[ -s ${ofile} ]]; then 
    efile_h="${ofile/.fni/-h.eps}"
    fni_hist -channel 0 -title "Z" -step 0.250 < ${ofile} > ${efile_h}
    if [[ -s ${efile_h} ]] ; then hfiles+=( ${efile_h} ); fi
  fi

  ofile=${out_slopes_fni};  if [[ -s ${ofile} ]]; then 
    efile_Xh="${ofile/.fni/X-h.eps}"
    fni_hist -channel 0 -title "dZ/DX" -step 0.250 < ${ofile} > ${efile_Xh}
    if [[ -s ${efile_Xh} ]] ; then hfiles+=( ${efile_Xh} ); fi

    efile_Yh="${ofile/.fni/Y-h.eps}"
    fni_hist -channel 1 -title "dZ/dY" -step 0.250 < ${ofile} > ${efile_Yh}
    if [[ -s ${efile_Yh} ]] ; then hfiles+=( ${efile_Yh} ); fi
  fi
fi

if [[ ${#pfiles[@]} -ne 0 ]]; then
  display -title '%f' -filter box -resize 'x480<' ${pfiles[@]}; status=$?
  if [[ ${status} -ne 0 ]]; then echo "** display exited with status = ${status} - aborted" 1>&2; exit 1; fi
fi

if [[ ${#efiles[@]} -ne 0 ]]; then
  evince ${efiles[@]}; status=$?
  if [[ ${status} -ne 0 ]]; then echo "** evince exited with status = ${status} - aborted" 1>&2; exit 1; fi
fi

if [[ ${#hfiles[@]} -ne 0 ]]; then
  evince ${hfiles[@]}; status=$?
  if [[ ${status} -ne 0 ]]; then echo "** evince exited with status = ${status} - aborted" 1>&2; exit 1; fi
fi
