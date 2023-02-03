#! /bin/bash
# Last edited on 2023-02-01 11:55:04 by stolfi

echo "=== do_regression_on_hist.sh =============================" 1>&2
echo "$@" 1>&2

prefix="$1"; shift      # File name minus the "-hdata.txt" tail.
basisName="$1"; shift   # Basis name ("DIFF", "LAPL", "HART", etc).
unitTerm="$1"; shift;   # True to include the unit term in the regression.
title="$1"; shift       # Plot title.

# Input files:
histDataFile="${prefix}-hdata.txt"    # File with quadratic terms binned by sharpness.
belNameFile="${prefix}-bnames.txt"    # File with basis element names.
termNameFile="${prefix}-tnames.txt"   # File with quadratic term names.

# Internal files:

regrDataFile="${prefix}-un${unitTerm}-rhdata.txt"  # Version of {histDataFile} modified for regression.

# Output files:

outFormFile="${prefix}-un${unitTerm}-hform.txt"  # Fitted formula for sharp as function of the binned quadratic terms.
outRegrFile="${prefix}-un${unitTerm}-hregr.txt"  # File with true and fitted sharpness.

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
cat ${histDataFile} \
  | gawk -v unitTerm=${unitTerm} \
      ' /^ *[0-9]/{ 
          ih = $1; sh = $2; 
          wt = sh;
          printf "%s %s %12.6f ", ih, sh, wt;
          for (kf = 3; kf <= NF; kf++) { printf " %s", $(kf) }
          if (unitTerm+0 > 0) { printf " 1"; }
          printf "\n";
        }
      ' \
  > ${regrDataFile}

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
 
