#! /bin/bash
# Last edited on 2023-04-28 18:51:37 by stolfi

sceneType="$1"; shift    # Scene type: "R" ramp only, "F" non-olap disks/balls, "T" olap.
NX="$1"; shift           # Image width.
NY="$1"; shift           # Image height.
patternDir="$1"; shift   # Directory where pattern images live, e.g "in/pgm-256x256".
pattern="$1"; shift      # Pattern name, e.g. "wavys-09".
zDep="$1"; shift         # Nominal depth of focus (pixels)
zFoc_min="$1"; shift     # Focus plane {Z} of first frame of stack.
zFoc_max="$1"; shift     # Max focus plane {Z} for stack frames.
zFoc_step="$1"; shift    # Increment in focus plane {Z} between stack frames.

# Input images:
#
#   in/pgm-256x256/${pattern}.pgm
# 
# Output images:
#
#   out/st{sceneType}-{size}-{pattern}/frame-sharp-c.pgm
#   out/st{sceneType}-{size}-{pattern}/frame-sharp-az.pgm
#   out/st{sceneType}-{size}-{pattern}/frame-sharp-dz.pgm
#   out/st{sceneType}-{size}-{pattern}/frame-sharp-sh.pgm (constant 1).
# 
#   out/st{sceneType}-{size}-{pattern}/frame-fd{zDep}-zf{zFoc}-cs.ppm
#   out/st{sceneType}-{size}-{pattern}/frame-fd{zDep}-zf{zFoc}-az.pgm
#   out/st{sceneType}-{size}-{pattern}/frame-fd{zDep}-zf{zFoc}-dz.pgm
#   out/st{sceneType}-{size}-{pattern}/frame-fd{zDep}-zf{zFoc}-sh.pgm
#
#   out/st{sceneType}-{size}-{pattern}/frame-sharp-show.ppm
#   out/st{sceneType}-{size}-{pattern}/frame-fd{zDep}-zf{zFoc}-show.ppm
#
# where {size} is "{NX}x{NY}" each formatted with "%04d", 
# and both {zDep} and {zFoc} are formatted with "%05.2f".

NP=5      # Number of pixel subsampling points in {X} and {Y}.
NR=40     # Min number of rays through each subsampling point.

size_FMT="`printf "%04dx%04d" ${NX} ${NY}`"
zDep_FMT="`printf "%05.2f" ${zDep}`"

# ----------------------------------------------------------------------
# Generate the images:

patternFile="${patternDir}/${pattern}.pgm"
outDir="out/st${sceneType}-${size_FMT}-${pattern}"
outPrefix="${outDir}/frame"
sharpImage="${outPrefix}-sharp-cs.ppm"

mkdir -pv ${outDir}
rm -f ${outPrefix}-sharp-*.{ppm,pgm}
rm -f ${outPrefix}-fd{zDep_FMT}-*.{ppm,pgm}

./test_mfok_make_image \
  -imgSize ${NX} ${NY} \
  -sceneType ${sceneType} \
  -pixSampling ${NP} \
  -dirSampling ${NR} \
  -zRange ${zFoc_min} ${zFoc_max} -zStep ${zFoc_step} \
  -focDepth ${zDep} \
  -patternFile ${patternFile} \
  -outPrefix ${outPrefix}
  
if [[ ! ( -s ${sharpImage} ) ]]; then echo "** images not generated" 1>&2; exit 1; fi

# ----------------------------------------------------------------------
# Display the images:

./show_mf_images.sh ${outPrefix} ${zDep_FMT}
