# Last edited on 2024-10-19 09:56:13 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_merge_images

JS_LIBS := \
  libimg.a \
  liblsq.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: check combine show

all: check

PATTERN := wavys-09
FOC_DEPTH := 6.0
SCENE_TYPE := F
NOISE := 0.04
FD_TAG := ${shell printf "df%08.4f" ${FOC_DEPTH}}

ZFOCS := -zRange 0.5 27.5 -zStep 3.0
ACTUAL_SHARP := F

IN_PREFIX := in/test-st${SCENE_TYPE}-${PATTERN}

OUT_PREFIX := out/test-st${SCENE_TYPE}-${PATTERN}-${FD_TAG}-ash${ACTUAL_SHARP}


check:  combine show

??MERGED_FILE := ${OUT_PREFIX}-cr.ppm

combine: ${MERGED_FILE}

${MERGED_FILE}: ${PROG} Makefile
	mkdir -p out
	rm -f ${OUT_PREFIX}*.ppm ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.txt
	./${PROG} \
            -inPrefix ${IN_PREFIX} \
            -focDepth ${FOC_DEPTH} \
            ${ZFOCS} \
            -noise ${NOISE} \
            -actualSharp ${ACTUAL_SHARP} \
            -outPrefix ${OUT_PREFIX}
            
show: ${MERGED_FILE}
	??display -title '%f' -filter Box -resize '400%' ${OUT_PREFIX}{-cr.ppm,-er.ppm,-az.pgm,-ez.pgm}
            
