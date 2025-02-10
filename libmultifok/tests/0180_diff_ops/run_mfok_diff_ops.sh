#! /bin/bash
# Last edited on 2025-01-30 10:30:22 by stolfi

echo "=== run_mfok_diff_ops.sh =============================" 1>&2
echo "$@" 1>&2

imageSet="$1"; shift;    # Name of image set file ("HARD", "SOSO", etc.).
basisType="$1"; shift;   # Name of local operato basis ("LAPL", "DIFF", etc.).
weightsType="$1"; shift; # Name of window weights distribution ("BINOMIAL", "GOLDEN", etc).
termSet="$1"; shift;     # Term subset to use ("ALL", etc.).
noise="$1"; shift;       # Noise level for local brightness/contrast normalization.

# ----------------------------------------------------------------------
# INPUT FILES
#
# For the given {termSet}, it reads a file 
#  
#   term-set/{basisType}-{termSet}.txt
#
# Each line should have the form
#
#   {kt} {wt[kt]} {termName[kt]}
#
# where an index {kt} in {0..NT-1}, {wt[kt]} is a non-negative weight
# (ignored) and {termname[kt]} is the formula of a quadratic term, a sum
# of product of pairs of basis elem names, like "DX*XX+DY*DY".
#
# For the given {imageSet}, this script reads a file 
#  
#   image-set/{imageSet}.txt
# 
# where each line has 
#
#   {sceneType} {size} {pattern} {zDep} {zFoc}

# identifying an input image.
# 
# For each specified combination of {sceneType},
# {size}, {pattern}, {zDep} and {zFoc} listed in the {imageSet} file,
# the program reads
#
#   {inPrefix}{frameTag}-cs.ppm
#   {inPrefix}{frameTag}-az.pgm
#   {inPrefix}{frameTag}-dz.pgm
#   {inPrefix}{frameTag}-sh.pgm
#
# where {inPrefix} is "in/st{sceneType}-{size}-{pattern}/frame" 
# and {frameTag} is "-zf{zFoc}"-df{zDep}
#
# OUTPUT FILES
#
# For the given parameters {basisType}, {weightsType}, {noise}, {imageSet}, the program writes
#
#   {outPrefix}-belnames.txt  Names of basis elements.
#   {outPrefix}-odata.txt     Pixel data incl. blurring indicator and coeffs before binning.
#   {outPrefix}-hdata.txt     Averages of {term[0..NT-1]}, binned by {shrp}.
#   {outPrefix}-hdata.png     Plot of terms, binned by sharpness{shrp}.
#
# where {outPrefix} is "out/bt{basisType}-wt${weightsType}-{termSet}-ns{noise}-ims{imageSet}/run"
# 
# And also for {unitTerm} 1 (add the "1" term to regression) and 0 (don't add that term):
#
#   {outUnPrefix}-oregr.txt    Output of linear regression, with {shrp} and fitted score per bin.
#   {outUnPrefix}-oregr.png    Plot of "{outUnPrefix}-oregr.txt".
#   {outUnPrefix}-oform.txt    Fitted regression formula to raw pixel data.
# 
#   {outUnPrefix}-hregr.txt    Output of linear regression on binned term averages, with given and fitted vals.
#   {outUnPrefix}-hregr.png    Plot of "{outUnPrefix}-hregr.txt".
#   {outUnPrefix}-hform.txt    Fitted regression formula to binned term averages.
#
# where {outUnPrefix} is "{outPrefix}-un{unitTerm}"
#
# Also for each given input combination of {sceneType},
# {size}, {pattern}, {zDep} and {zFoc} included in the {imageSet}:
#
#   {outImagePrefix}{frameTag}-av.pgm   Window average of grayed "-cs.ppm" input image.
#   {outImagePrefix}{frameTag}-gv.pgm   Window gradient of grayed "-cs.ppm" input image.
#   {outImagePrefix}{frameTag}-dv.pgm   Window deviation of grayed "-cs.ppm" input image.
#   {outImagePrefix}{frameTag}-nr.pgm   Locally normalized version of grayed "-cs.ppm" input image.
#   {outImagePrefix}{frameTag}-mk.pgm   Mask showing which pixels were considered in the binned averages.
#
# where {outImagePrefix} is "{outPrefix}-img-st{sceneType}-{size}-{pattern}".
#
# Also for given {outFrameTag} and each basis elem index {kb}
#
#   {outImagePrefix}{frameTag}-{kb}-bq.pgm   Value of {coeff[kb]^2} at each pixel.
#
# Also for given {outFrameTag} and each quadratic term index {kt}
#
#   {outImagePrefix}{frameTag}-{kt}-tm.pgm   Value of {term[kt]} at each pixel.
#
# Finally it writes a single file
#
#   out/bt${basisType}-wt${weightsType}-ns${noise_TAG}-alltcoefs.txt
#
# which has the term coefficients of all linear regressions, for all {imgSet} and {unitTerm},
# side by side.

noise_TAG="`printf "%05.2f" ${noise}`"

inFolder="in"

formFiles=( ) # Files with term coeffcients fitted to binned or pixel data.

# Get the list of images to analyze:
wc -l image-set/${imageSet}.txt
imageOptions=( `cat image-set/${imageSet}.txt | sed -e 's:[#].*$::g' -e '/[A-Z]/s:^:-image :g'` )
no=${#imageOptions[@]}
echo "no = ${no}" 1>&2
ko=0
while [[ ${ko} -lt ${no} ]]; do
  echo "[${imageOptions[@]:${ko}:7}]" 1>&2
  ko=$(( ${ko} + 7 ))
done

outFolder="out/bt${basisType}-wt${weightsType}-${termSet}-ns${noise_TAG}-ims${imageSet}"
outPrefix="${outFolder}/run"

# Get the list of quadratic terms descriptions:
termSetFile="term-set/${basisType}-${termSet}.txt"

mkdir -pv ${outFolder}
rm -fv ${outPrefix}*.ppm ${outPrefix}*.pgm ${outPrefix}*.txt
./test_mfok_diff_ops \
  -inFolder ${inFolder} \
  ${imageOptions[@]} \
  -basisType ${basisType} \
  -weightsType ${weightsType} \
  -termSetFile ${termSetFile} \
  -noise ${noise} \
  -outPrefix ${outPrefix}

# Get the basis element names and count:
belNameFile="${outPrefix}-bnames.txt"      # File with basis element names.
if [[ ! ( -s ${belNameFile} ) ]]; then echo "** missing file \"${belNameFile}\"" 1>&2; exit 1; fi
belName=( `cat ${belNameFile}` )
nb=${#belName[@]}
echo "found ${nb} basis elements" 1>&2 
if [[ ${nb} -gt 9 ]]; then echo "** too many coeffs" 1>&2; exit 1; fi

# Get the quadratic term names and count:
termNameFile="${outPrefix}-tnames.txt"
if [[ ! ( -s ${termNameFile} ) ]]; then echo "** missing file \"${termNameFile}\"" 1>&2; exit 1; fi
termName=( `cat ${termNameFile} | cut -f 3` )
nt=${#termName[@]}
echo "found ${nt} quadratic elements" 1>&2 
if [[ ${nt} -gt 45 ]]; then echo "** too many coeffs" 1>&2; exit 1; fi

# File with one line per pixel:
pixDataFile="${outPrefix}-odata.txt"
if [[ -s ${pixDataFile} ]]; then

  wc -l ${pixDataFile} 1>&2

  histDataFile="${outPrefix}-hdata.txt"      # File with histogram data for plot and regression.
  echo "condensing quadratic terms by {shrp} bin -> ${histDataFile} ..." 1>&2 
  cat ${pixDataFile} \
    | ./average_terms_by_hrad.gawk \
        -v nb=${nb} \
        -v nt=${nt} \
    > ${histDataFile}
  wc -l ${histDataFile}

  coeffPlotTitle="Image set ${imageSet}  term set ${basisType}-${termSet}  noise ${noise}"
  echo "plot the averages of quadratic terms binned by {shrp} -> ${coeffPlotTitle} ..." 1>&2 
  ./plot_terms_by_hrad2.sh ${outPrefix} ${basisType} ${termSet} "${coeffPlotTitle}"

  echo "showing the quadratic terms term images for the firsts input image ..." 1>&2 
  ./show_term_images.sh ${inFolder} ${outPrefix} ${imageOptions[@]:1:6}

  # Do regression on {shrp} as func of bin-averaged quadratic terms with and without the "1" term:
  for unitTerm in 0 1 ; do

    if [[ ${unitTerm} -gt 0 ]]; then unitTerm_TITLE="+UNIT"; else unitTerm_TITLE=""; fi
    regrPlotTitle="Image set ${imageSet} terms ${basisType}-${termSet}${unitTerm_TITLE} noise ${noise}"

    histFormFile=${outPrefix}-un${unitTerm}-hform.txt
    echo "performing regression on binned data -> ${histFormFile} ..." 1>&2 
    do_regression_on_binned_data.sh ${outPrefix} ${basisType} ${unitTerm} "${regrPlotTitle} (BINNED)"
    if [[ ! (-s ${histFormFile} ) ]]; then
      echo "** binned data regression failed - no file \"${histFormFile}\"" 1>&2 ; exit 1
    fi
    formfiles+=( ${histFormFile} )

    # Regression on raw pixel data:

    pixFormFile=${outPrefix}-un${unitTerm}-oform.txt
    if [[ ${unitTerm} -gt 0 ]]; then unitTerm_TITLE="+UNIT"; else unitTerm_TITLE=""; fi
    echo "performing regression on raw pixel data -> ${pixFormFile} ..." 1>&2 
    do_regression_on_pix_data.sh ${outPrefix} ${basisType} ${unitTerm} "${regrPlotTitle} (PIXEL)"
    if [[ ! (-s ${pixFormFile} ) ]]; then
      echo "** pixel data regression failed - no file \"${pixFormFile}\"" 1>&2 ; exit 1
    fi
    formfiles+=( ${pixFormFile} )

  done # Loop on ${unitTerm}
else
  echo "** mf_0180_diff_ops failed - no file \"${pixDataFile}\"" 1>&2 ; exit 1
fi

jointTermsFile="${outPrefix}-alltcoefs.txt"
paste_formulas.gawk "${formfiles[@]}" "${Termfiles[@]}" > ${jointTermsFile}
cat ${jointTermsFile}
