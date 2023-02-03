#! /bin/bash
# Last edited on 2023-02-01 19:00:19 by stolfi

echo "=== run_0180_diff_ops.sh =============================" 1>&2
echo "$@" 1>&2

basisName="$1"; shift;   # Name of local operato basis ("LAPL", "DIFF", etc.)
termSet="$1"; shift;     # Term subset to use ("ALL", etc).
noise="$1"; shift;       # Noise level for local brightness/contrast normalization.

# ----------------------------------------------------------------------
# INPUT FILES
#
# For the given {termSet}, it reads a file 
#  
#   term-set/{basisName}-{termSet}.txt
#
# Where each line has an index {kt} in {0..NT-1} and the formula
# of a quadratic term, a sum of product of pairs of basis elem names, like
# "DX*XX+DY*DY".
#
# For each {imageSet} specified below, this script reads a file 
#  
#   image-set/{imageSet}.txt
# 
# where each line has {sceneType} {size}, {pattern}, {zDep} and {zFoc}
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
# and {frameTag} is "-fd{zDep}-zf{zFoc}"
#
# OUTPUT FILES
#
# For the given parameters {basisName}, {noise}, {imageSet}, the program writes
#
#   {outPrefix}-belnames.txt  Names of basis elements.
#   {outPrefix}-odata.txt     Pixel data incl. sharpness and coeffs before binning.
#   {outPrefix}-hdata.txt     Averages of squared {term[0..NT-1]}, binned by sharpness.
#   {outPrefix}-hdata.png     Plot of squared coeffs, binned by sharpness.
#
# where {outPrefix} is "out/bt{basisName}-{termSet}-ns{noise}-ims{imageSet}/run"
# 
# And also for {unitTerm} 1 (add the "1" term to regression) and 0 (don't add that term):
#
#   {outUnPrefix}-ordata.txt   Variant of {outPrefix}-odata.txt} possibly with extra unit term colum.
#   {outUnPrefix}-oregr.txt    Output of linear regression, with sharpness and fitted score per bin.
#   {outUnPrefix}-oregr.png    Plot of "{outUnPrefix}-oregr.txt".
#   {outUnPrefix}-oform.txt    Fitted regression formula to raw pixel data.
# 
#   {outUnPrefix}-hrdata.txt   Variant of {outPrefix}-hdata.txt} possibly with extra unit term colum.
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
#   out/bt${basisName}-ns${noise_TAG}-alltcoefs.txt
#
# which has the term coefficients of all linear regressions, for all {imgSet} and {unitTerm},
# side by side.

noise_TAG="`printf "%05.2f" ${noise}`"

inDir="in"

formFiles=( ) # Files with term coeffcients fitted to binned or pixel data.

# imageSets=( BEST SOSO HARD )
imageSets=( HARD )

for imageSet in ${imageSets[@]} ; do

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

  outDir="out/bt${basisName}-${termSet}-ns${noise_TAG}-ims${imageSet}"
  outPrefix="${outDir}/run"

  # Get the list of quadratic terms descriptions:
  termSetFile="term-set/${basisName}-${termSet}.txt"
  if [[ ! ( -s ${termSetFile} ) ]]; then echo "** missing file \"${termSetFile}\"" 1>&2; exit 1; fi
  termOptions=( `cat ${termSetFile} | gawk '/^ *[0-9]/ { print "-term", $2; }'` )
  echo "termOptions = ${termOptions[*]}" 1>&2

  mkdir -pv ${outDir}
  rm -fv ${outPrefix}*.ppm ${outPrefix}*.pgm ${outPrefix}*.txt
  ./mf_0180_diff_ops \
    -inDir ${inDir} \
    ${imageOptions[@]} \
    -basisName ${basisName} \
    ${termOptions[@]} \
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
  termName=( `cat ${termNameFile} | cut -f 2` )
  nt=${#termName[@]}
  echo "found ${nt} quadratic elements" 1>&2 
  if [[ ${nt} -gt 45 ]]; then echo "** too many coeffs" 1>&2; exit 1; fi

  # File with one line per pixel:
  pixDataFile="${outPrefix}-odata.txt"
  if [[ -s ${pixDataFile} ]]; then

    wc -l ${pixDataFile} 1>&2

    # Condense the quadratic terms by sharpness bin:
    histDataFile="${outPrefix}-hdata.txt"      # File with histogram data for plot and regression.

    cat ${pixDataFile} \
      | ./filter_pixel_data.gawk \
      | ./average_terms_by_sharp.gawk \
          -v nb=${nb} \
          -v nt=${nt} \
      > ${histDataFile}
    wc -l ${histDataFile}

    # Plot the averages of quadratic terms binned by sharpness:
    coeffPlotTitle="Image set ${imageSet}  term set ${basisName}-${termSet}  noise ${noise}"
    ./plot_terms_by_sharp.sh ${outPrefix} ${basisName} ${termSet} "${coeffPlotTitle}"

    # Show images with the quadratic terms for the first image: 
     ./show_term_images.sh ${inDir} ${outPrefix} ${imageOptions[@]:1:6}

    # Do regression on sharpness as func of bin-averaged coeffs squared with and without the "1" term:
    for unitTerm in 0 1 ; do
    
      if [[ ${unitTerm} -gt 0 ]]; then unitTerm_TITLE="+UNIT"; else unitTerm_TITLE=""; fi
      regrPlotTitle="Image set ${imageSet} terms ${basisName}-${termSet}${unitTerm_TITLE} noise ${noise}"
      
      # Regression on binned data:
    
      # File with formula fitted to histogram:
      histFormFile=${outPrefix}-un${unitTerm}-hform.txt

      do_regression_on_hist.sh ${outPrefix} ${basisName} ${unitTerm} "${regrPlotTitle} (BINNED)"
      if [[ ! (-s ${histFormFile} ) ]]; then
        echo "** binned data regression failed - no file \"${histFormFile}\"" 1>&2 ; exit 1
      fi
      formfiles+=( ${histFormFile} )
    
      # Regression on raw pixel data:
    
      # File with formula fitted to pixel data:
      pixFormFile=${outPrefix}-un${unitTerm}-oform.txt
      if [[ ${unitTerm} -gt 0 ]]; then unitTerm_TITLE="+UNIT"; else unitTerm_TITLE=""; fi
      do_regression_on_pix_data.sh ${outPrefix} ${basisName} ${unitTerm} "${regrPlotTitle} (PIXEL)"
      if [[ ! (-s ${pixFormFile} ) ]]; then
        echo "** pixel data regression failed - no file \"${pixFormFile}\"" 1>&2 ; exit 1
      fi
      formfiles+=( ${pixFormFile} )

    done # Loop on ${unitTerm}
  else
    echo "** mf_0180_diff_ops failed - no file \"${pixDataFile}\"" 1>&2 ; exit 1
  fi

done # Loop on ${imageSet}

jointTermsFile="out/bt${basisName}-ns${noise_TAG}-alltcoefs.txt"
paste_formulas.gawk "${formfiles[@]}" "${Termfiles[@]}" > ${jointTermsFile}
