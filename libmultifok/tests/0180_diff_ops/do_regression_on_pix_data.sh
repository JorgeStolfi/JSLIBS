#! /bin/bash
# Last edited on 2023-02-01 18:10:30 by stolfi

echo "=== ${0/*\/} =============================" 1>&2
echo "$@" 1>&2

prefix="$1"; shift      # File name minus the "-odata.txt" tail.
basisName="$1"; shift   # Basis name ("DIFF", "LAPL", "HART", etc).
unitTerm="$1"; shift;   # True to include the unit term in the regression.
title="$1"; shift       # Plot title.

# Input files:
pixDataFile="${prefix}-odata.txt"     # File with per-pixel sharpness and term data.
belNameFile="${prefix}-bnames.txt"    # File with basis element names.
termNameFile="${prefix}-tnames.txt"   # File with quadratic term names.

# Internal files:

regrDataFile="${prefix}-un${unitTerm}-rodata.txt"  # Version of {pixDataFile} modified for regression.

# Output files:

outFormFile="${prefix}-un${unitTerm}-oform.txt"  # Fitted formula for sharp as function of the quadratic terms.
outRegrFile="${prefix}-un${unitTerm}-oregr.txt"  # File with true and fitted sharpness.

# Get the basis element names:
belName=( `cat ${belNameFile}` )
nb=${#belName[@]}
echo "found ${nb} basis elements" 1>&2 
if [[ ${nb} -gt 9 ]]; then echo "** too many coeffs" 1>&2; exit 1; fi

# Get the quadratic term names:
termName=( `cat ${termNameFile}` )
if [[ ${unitTerm} -ne 0 ]]; then
  termName+=( "1" )
  unitX=" (including unit)"
fi
nt=${#termName[@]}
echo "found ${nt} quadratic terms${unitX}" 1>&2 
echo "termName = ${termName[*]}" 1>&2 

# Prepare the data for regression:
echo "generating the regression input file ${regrDataFile}..." 1>&2
cat ${pixDataFile} \
  | filter_pixel_data.gawk \
  | gawk \
      -v unitTerm=${unitTerm} \
      -v nb=${nb} \
      -v nt=${nt} \
      ' /^ *P[0-9]/{ 
          if (NF + unitTerm  != 6 + nb + nt) { printf "BUG\n" > "/dev/stderr"; exit(1); }
          pi = $1; va = $2; vd = $3; sh = $4; za = $5; zd = $6; 
          # wp = sh;
          wp = 1.0;
          printf "%s %s %12.6f ", pi, sh, wp;
          for (kf = 7+nb; kf <= NF; kf++) { printf " %s", $(kf) }
          if (unitTerm+0 > 0) { printf " 1"; }
          printf "\n";
        }
      ' \
  > ${regrDataFile}
  
if [[ ! ( -s ${regrDataFile} ) ]]; then
  echo "** failed to create the regression data file" 1>&2; exit 1
fi

# Do a regression on the modified histogram file:
linear_fit \
    -terms ${nt} \
    -weighted T \
    -termNames "${termName[@]}" \
    -writeFormula ${outFormFile} \
  < ${regrDataFile} \
  > ${outRegrFile}

# Plot summary of regression:
plot_regression_result.sh SHOW ${outRegrFile/.txt/} "${title}" "Sharpness" "Fitted"
 
