#! /bin/bash
# Last edited on 2023-01-25 07:54:04 by stolfi

prefix="$1"; shift      # File name minus the "-terma.txt" and "-rdata.txt" tail.

tmp="/tmp/$$"

termsFile="${prefix}-terms.txt" # File with term names and weights.
dataFile="${prefix}-rdata.txt"  # Raw data file.

termNames=( `cat ${termsFile} | sed -e 's:[#].*$::g' | gawk '//{ print $1; }'` )
nt=${#termNames[@]}
echo "expecting ${nt} terms in regression data file" 1>&2

# Assumes that each line of ${dataFile} corresponds to a pixel and has fields
#   {pixelID} {sharp^exp} {wp} {term[0]} ... {term[NT-1]}
# where {wp} is the weight of the pixel, {sharp} is the actual image sharpness 
# in the pixel, {exp} is 0 or 1, and {term[0..NT-1]} are the terms (local 
# operator values) that are to be used to estimate the {sharp} value.

outRegrFile="${prefix}-regr.txt"   # Regression result file.
outFormFile="${prefix}-form.txt"   # Fitted formula file.

# Do the regression: */

linear_fit \
    -terms ${nt} \
    -weighted T \
    -termNames "${termNames[@]}" \
    -writeFormula ${outFormFile} \
  < ${dataFile} \
  > ${outRegrFile}
  
