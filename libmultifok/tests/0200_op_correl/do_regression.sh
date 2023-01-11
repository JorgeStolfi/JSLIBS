#! /bin/bash
# Last edited on 2023-01-09 03:39:58 by stolfi

fname="$1"; shift   # File name minus the "-b{B}-vals.txt" tail.
btype="$1"; shift   # "N", "D", or "H".
nw="$1"; shift      # Window size.
sqrOnly="$1"; shift # 1 to generate only the squared terms, 0 all cross products too.

if [[ ${nw} -ne 3 ]]; then echo "** window size must be 3" 1>&2; exit 1; fi
nb=$(( $nw * $nw ))

tmp="/tmp/$$"

inNamesFile="${fname}-b${btype}-names.txt" # File with basis element names.
inCoeffFile="${fname}-b${btype}-vals.txt"  # Raw data file.

tmpTermsFile="${tmp}-tdata.txt"             # Temporary data file with quadratic terms.
tmpNamesFile="${tmp}-tnames.txt"            # Temporary data file with names of quadratic terms.
outRegrFile="${fname}-b${btype}-regr.txt"   # Regression result file.
outFormFile="${fname}-b${btype}-form.txt"   # Fitted formula file.

# Get the number of basis elements and their names:

names=( `cat ${inNamesFile}` )
nb=${#names[@]}

echo "found ${nb} basis elements = ${names[@]}" 1>&2 

# Generating the term names file:

gawk \
    -v names="${names[*]}" \
    -v sqrOnly=${sqrOnly} \
    ' BEGIN {
        nb = split(names, bnam);
        printf "generating the terms found %d basis elements\n", nb > "/dev/stderr"
        for (ib = 1; ib <= nb; ib++) {
          dblim = (sqrOnly+0 != 0 ? 1 : ib);
          for (db = 0; db < dblim; db++) {
            jb = ib-db;
            tn = ( bnam[ib] "*" bnam[jb] )
            printf "%s\n", tn;
          }
        }
        printf "1\n";
      }
    ' \
  > ${tmpNamesFile}
  
if [[ ! ( -s ${tmpNamesFile} ) ]]; then exit 1; fi
  
if [[ ${sqrOnly} -ne 0 ]]; then
  nt=$(( ${nb} + 1 ))
else
  nt=$(( ${nb} * ( ${nb} + 1) / 2 + 1 )) 
fi
echo "should have ${nt} terms; checking" 1>&2
wc -l ${tmpNamesFile} 1>&2

# Generate the temporary data file with the terms values:

cat ${inCoeffFile} \
  | gawk \
      -v nb="${nb}" \
      -v sqrOnly=${sqrOnly} \
      ' // { 
          printf "P%03d.%03d %s ", $1, $2, $4; 
          for (ib = 1; ib <= nb; ib++) {
            dblim = (sqrOnly+0 != 0 ? 1 : ib);
            for (db = 0; db < dblim; db++) {
              jb = ib-db;
              tv = $(4+ib) * $(4+jb);
              printf " %16.12f", tv;
            }
          }
          printf " 1\n";
        }
      ' \
  > ${tmpTermsFile}
  
if [[ ! ( -s ${tmpNamesFile} ) ]]; then exit 1; fi

# Do the quadratic regression:

linear_fit \
    -terms ${nt} \
    -termNames `cat ${tmpNamesFile}` \
    -writeFormula ${outFormFile} \
  < ${tmpTermsFile} \
  > ${outRegrFile}
  
