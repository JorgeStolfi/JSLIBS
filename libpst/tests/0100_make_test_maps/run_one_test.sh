#! /bin/bash
# Last edited on 2025-01-25 17:13:42 by stolfi

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

outDir=out
sizeTag="`printf "%04dx%04d" ${nx} ${ny}`"
outPrefix="${outDir}/${func_num}-${func_name}-${noisy}-${sizeTag}"

rm -f ${outPrefix}-*.{fni,ppm,pgm,eps}
set -x
${PROGDIR}/${PROG} \
  -function ${func_num} \
  -size ${nx} ${ny} \
  -numGrad Y \
  -noiseG ${noiseG} \
  -outDir ${outDir}
set +x

ls -d ${outPrefix}-*.fni

out_heights_fni="${outPrefix}-Z.fni"
out_slopes_n_fni="${outPrefix}-Gn.fni"
out_slopes_a_fni="${outPrefix}-Ga.fni"
out_slopes_d_fni="${outPrefix}-Gd.fni"
out_normals_fni="${outPrefix}-N.fni"

view_fni=1
make_pnm=1
make_hist=0
make_plot=0

# ======================================================================
# Viewing maps with fni_view

function view(){
  ofile="$1"; shift;
  txSelf="$1"; shift
  
  scale=`echo "sqrt(${nx}^2+${ny}^2)/2" | bc -lq`
  scale=`printf "%.2f" "${scale}"`
  if [[ ${txSelf} -ne 0 ]]; then txOps=( -txChannels 0 1 2 ); else txOps=(); fi
  if [[ -s ${ofile} ]]; then
    set -x
    fni_view -channel 0 -hist T -range auto -scale ${scale} ${txOps[@]} -verbose ${ofile}
    set +x
  else
    echo "** did not find file ${ofile}" 1>&2
  fi
}

if [[ ${view_fni} -gt 0 ]]; then
  view ${out_heights_fni}   0
  view ${out_slopes_a_fni}  0
  view ${out_slopes_n_fni}  0
  view ${out_slopes_d_fni}  0
  view ${out_normals_fni}   1
fi

# ======================================================================
# Converting maps to ".ppm" and/or ".pgm"

pfiles=()

function topnm(){
  ofile="$1"; shift
  chans="$1"; shift
  tag="$1"; shift
  
  chans=(`echo ${chans} | tr ',' ' '`)
  
  if [[ ${#chans[@]} -eq 1 ]]; then
    ext="pgm"; chop=( -channel ${chans[0]} )
  elif [[ ${#chans[@]} -eq 3 ]]; then
    ext="ppm"; chop=( -channels ${chans[@]} )
  else
    echo "** invalid {chans}" 1>&2; exit 1;
  fi
  if [[ -s ${ofile} ]]; then 
    pfile="${ofile/.fni/${tag}.${ext}}"
    fni_to_pnm ${chop[@]} -center 0 -uniform < ${ofile} > ${pfile}
    if [[ -s ${pfile} ]]; then pfiles+=( ${pfile} ); fi
  else
    echo "** file ${ofile} not found" 1>&2
  fi
}

if [[ ${make_pnm} -gt 0 ]]; then
  
  topnm ${out_heights_fni} 0    ''
  topnm ${out_heights_fni} 1    'W'
  
  topnm ${out_slopes_a_fni} 0,1,0 ''
  topnm ${out_slopes_a_fni} 2     'W'
  
  topnm ${out_slopes_n_fni} 0,1,0 ''
  topnm ${out_slopes_n_fni} 2     'W'
  
  topnm ${out_slopes_d_fni} 0,1,0 ''
  topnm ${out_slopes_d_fni} 2     'W'
  
  topnm ${out_normals_fni} 0,1,2 ''
  topnm ${out_normals_fni} 3     'W'

  if [[ ${#pfiles[@]} -ne 0 ]]; then
    display -title '%f' -filter box -resize 'x480<' ${pfiles[@]}; status=$?
    if [[ ${status} -ne 0 ]]; then echo "** display exited with status = ${status} - aborted" 1>&2; exit 1; fi
  fi
fi  

# ======================================================================
# Converting maps to EPS plots:

efiles=()

function plotit(){
  ofile="$1"; shift
  chan="$1"; shift
  tag="$1"; shift
  title="$1"; shift
  
  if [[ -s ${ofile} ]]; then 
    efile="${ofile/.fni/${tag}.eps}"
    fni_plot.sh -channel ${chan} -title "${title}" < ${ofile} > ${efile}
    if [[ -s ${efile} ]]; then efiles+=( ${efile} ); fi
  else
    echo "** file ${ofile} not found" 1>&2
  fi
}

if [[ ${make_plot} -gt 0 ]]; then

  plotit ${out_heights_fni} 0 '' "Z"

  plotit ${out_slopes_a_fni} 0 'X' "dZ/dX"
  plotit ${out_slopes_a_fni} 1 'Y' "dZ/dY"

  plotit ${out_slopes_n_fni} 0 'X' "dZ/dX"
  plotit ${out_slopes_n_fni} 1 'Y' "dZ/dY"

  plotit ${out_slopes_d_fni} 0 'X' "dZ/dX"
  plotit ${out_slopes_d_fni} 1 'Y' "dZ/dY"

  plotit ${out_normals_fni} 0 'X' "nrm.X"
  plotit ${out_normals_fni} 0 'Y' "nrm.Y"
  plotit ${out_normals_fni} 0 'Z' "nrm.Z"

  if [[ ${#efiles[@]} -ne 0 ]]; then
    evince ${efiles[@]}; status=$?
    if [[ ${status} -ne 0 ]]; then echo "** evince exited with status = ${status} - aborted" 1>&2; exit 1; fi
  fi

fi

# ======================================================================
# Computing map histograms:

hfiles=()

function tohist(){
  ofile="$1"; shift
  chan="$1"; shift
  tag="$1"; shift
  title="$1"; shift
  
  if [[ -s ${ofile} ]]; then 
    hfile="${ofile/.fni/${tag}-hist.eps}"
    fni_hist -channel ${chan} -step 0.250 -title "${title}" < ${ofile} > ${hfile}
    if [[ -s ${hfile} ]]; then hfiles+=( ${hfile} ); fi
  else
    echo "** file ${ofile} not found" 1>&2
  fi
}

if [[ ${make_hist} -gt 0 ]]; then

  tohist ${out_heights_fni} 0 '' "Z"

  tohist ${out_slopes_a_fni} 0 'X' "dZ/dX"
  tohist ${out_slopes_a_fni} 1 'Y' "dZ/dY"

  tohist ${out_slopes_n_fni} 0 'X' "dZ/dX"
  tohist ${out_slopes_n_fni} 1 'Y' "dZ/dY"

  tohist ${out_slopes_d_fni} 0 'X' "dZ/dX"
  tohist ${out_slopes_d_fni} 1 'Y' "dZ/dY"

  if [[ ${#hfiles[@]} -ne 0 ]]; then
    evince ${hfiles[@]}; status=$?
    if [[ ${status} -ne 0 ]]; then echo "** evince exited with status = ${status} - aborted" 1>&2; exit 1; fi
  fi
fi
