#! /bin/bash
# Last edited on 2023-04-28 11:32:54 by stolfi

echo "=== do_regression_on_binned_data.sh =============================" 1>&2
echo "$@" 1>&2

prefix="$1"; shift      # File name minus the "-hdata.txt" tail.
basisName="$1"; shift   # Basis name ("DIFF", "LAPL", "HART", etc).
unitTerm="$1"; shift;   # True to include the unit term in the regression.
title="$1"; shift       # Plot title.

# Input files:
histDataFile="${prefix}-hdata.txt"    # File with quadratic terms binned by sharpness.
belNameFile="${prefix}-bnames.txt"    # File with basis element names.
termNameFile="${prefix}-tnames.txt"   # File with quadratic term names.

# Output files:

outFormFile="${prefix}-un${unitTerm}-hform.txt"  # Fitted formula for sharp as function of the binned quadratic terms.
outRegrFile="${prefix}-un${unitTerm}-hregr.txt"  # File with true and fitted sharpness.

# Internal files:

regrDataFile="${prefix}-un${unitTerm}-rtmhist.txt"  # Version of {histDataFile} modified for regression.

# Get the basis element names:
belName=( `cat ${belNameFile}` )
nb=${#belName[@]}
echo "found ${nb} basis elements" 1>&2 
if [[ ${nb} -gt 9 ]]; then echo "** too many coeffs" 1>&2; exit 1; fi

# Get the quadratic term names:
termName=( `cat ${termNameFile}` )
nt=${#termName[@]}
echo "found ${nt} quadratic terms" 1>&2 
echo "termName = ${termName[*]}" 1>&2 

# Preparing data file for regression:
echo "converting ${histDataFile} to regression data file ${regrDataFile}..." 1>&2
cat ${histDataFile} \
  | gawk  \
      ' /^ *[0-9]/{ 
          pixid = $1; hrad = $2; 
          wt = 1.0/hrad;   # Weight equal to "true" sharpness.
          hrad2 = hrad*hrad;
          printf "%s %12.6f %12.6f ", pixid, hrad2, wt;
          for (kf = 3; kf <= NF; kf++) { printf " %14.8f", $(kf) }
          printf "\n";
        }
      ' \
  > ${regrDataFile}

# Do the regression of {hrad^2} against the quadtratic terms:
echo "performing regression on ${regrDataFile}..." 1>&2
linear_fit \
    -terms ${nt} \
    -weighted T \
    -unitTerm ${unitTerm} \
    -termNames "${termName[@]}" \
    -verbose T \
    -writeFormula ${outFormFile} \
  < ${regrDataFile} \
  > ${outRegrFile}

# Plot summary of regression:
plot_regression_result.sh SHOW ${outRegrFile/.txt/} "${title}" "hrad2" "Fitted"
 
