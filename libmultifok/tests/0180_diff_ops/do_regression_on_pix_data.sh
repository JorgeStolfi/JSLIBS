#! /bin/bash
# Last edited on 2024-10-16 13:47:56 by stolfi

echo "=== ${0/*\/} =============================" 1>&2
echo "$@" 1>&2

prefix="$1"; shift      # File name minus the "-odata.txt" tail.
basisType="$1"; shift   # Basis name ("DIFF", "LAPL", "HART", etc).
unitTerm="$1"; shift;   # True to include the unit term in the regression.
title="$1"; shift       # Plot title.

# Input files:
pixDataFile="${prefix}-odata.txt"     # File with per-pixel blurring indicator and term data.
belNameFile="${prefix}-bnames.txt"    # File with basis element names.
termNameFile="${prefix}-tnames.txt"   # File with quadratic term names.

# Output files:

outFormFile="${prefix}-un${unitTerm}-oform.txt"  # Fitted formula for {shrp} as function of the quadratic terms.
outRegrFile="${prefix}-un${unitTerm}-oregr.txt"  # File with true and fitted {shrp} indicator.

# Internal files:

regrDataFile="${prefix}-un${unitTerm}-rtmpix.txt"  # Version of {pixDataFile} modified for regression.

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
echo "preparing regression data file ${regrDataFile}..." 1>&2

cat ${pixDataFile} \
  | gawk  \
      -v nb=${nb} \
      -v nt=${nt} \
      ' /^ *P[0-9]/{ 
          if (NF != 7 + nb + nt) { printf "BUG\n" > "/dev/stderr"; exit(1); }
          pixid = $1; vavg = $2; vgrd = $3; vdev = $4; shrp = $5; zrav = $6; zdev = $7; 
          wt = sharp??;   # Weight equal to "true" sharpness.
          # wt = 1.0;  # Uniform weight.
          hrad = (1.0/sharp??);
          hrad2 = hrad*hrad;
          printf "%s %12.6f %12.6f ", pixid, hrad2, wt;
          for (kt = 0; kt < nt; kt++) { kf = 8 + nb + kt; printf " %s", $(kf) }
          printf "\n";
        }
      ' \
  | sort -b -k2g \
  > ${regrDataFile}

# Do the regression:
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
 
