# Last edited on 2024-10-28 15:14:00 by stolfi

SHELL := /bin/bash

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_scene_make_frame

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: \
  all-stacks \
  all-stacks-show-images \
  all-stacks-plot-profiles \
  single-stack \
  single-stack-show-images \
  single-stack-plot-profiles \
  single-stack-clean

all: \
  all-stacks \
  all-stacks-show-images \
  all-stacks-plot-profiles 

# Image dimensions:
# IMAGE_NX := 640
# IMAGE_NY := 480

# IMAGE_NX := 160
# IMAGE_NY := 120

IMAGE_NX := 320
IMAGE_NY := 240

# Subpixel and ray sampling:
PIX_SAMPLING := 2  # {HS} in code.  There will be {(2*HS+1)^2} sampling points per pixel.
DIR_SAMPLING := 3  # {KR} in code.  There will be at least {KR} rays per sampling point.

# {TEXTURE}:{MULTIPLICITY}:{SCENE_TYPE}
# where {MULTIPLICITY} is "S" for single frame, "M" for frame stack,
# and {SCENE_TYPE} is "R" for ramp only, "F" for non-overlapping objs,
# "T" for overlapping.

# Quick test to debug the program (single frame of ramp):
QUICK_TESTS := \
  noise01:S:R
  
# Basic test for multifocus stereo (ramp with fine texture):
BASIC_TESTS := \
  noise01:M:R
  
OBJECT_TESTS := \
  noise01:M:F

SOME_TESTS := \
  cloth84:S:R \
  fiber01:S:R \
  fiber01:S:R \
  grave02:S:R  \
  grittie:S:R\
  noise01:S:R\
  pebbl87:S:R \
   \
  cloth84:S:F \
  grave02:S:F \
  grittie:S:F \
  noise01:S:F \
   \
  fiber01:M:F \
  grittie:M:F \
  noise01:M:F \

MORE_TESTS := \
  grave02:S:R \
  grittie:S:R \
   \
  decor82:S:F \
  falls87:S:F \
  fiber01:S:F \
  gecko02:S:F \
  gravl87:S:F \
  husks01:S:F \
  india81:S:F \
  leafy09:S:F \
  mosic96:S:F \
  pebbl87:S:F \
  rebar03:S:F \
  sqirt01:S:F \
  tampi86:S:F \
  tarmc81:S:F \
  trees07:S:F \
  wavys33:S:F \
  wheat81:S:F \


ALL_TESTS := \
  ${SOME_TESTS} \
  ${MORE_TESTS}

# TESTS := ${QUICK_TESTS}
# TESTS := ${BASIC_TESTS}
TESTS := ${OBJECT_TESTS}

all-stacks: 
	for targs in ${TESTS} ; do \
	  ${MAKE} TARGS=$${targs} single-stack ; \
        done
        
all-stacks-show-images: 
	for targs in ${TESTS} ; do \
	  ${MAKE} TARGS=$${targs} single-stack-show-images ; \
        done

all-stacks-plot-profiles: 
	for targs in ${TESTS} ; do \
	  ${MAKE} TARGS=$${targs} single-stack-plot-profiles ; \
        done

######################################################################
# recursive section - caller must define ${TARGS}
ifneq "/${TARGS}" "/"

TEST_ARGS :=  ${subst :, , ${TARGS}}
PATTERN :=    ${word 1, ${TEST_ARGS}}
ZOPT :=       ${word 2, ${TEST_ARGS}}
SCENE_TYPE := ${word 3, ${TEST_ARGS}}

Z_DEP := 12.0
ZFOC_STEP := ${shell v=`echo '${Z_DEP}/2.0' | bc -lq` ; printf "%.3f" $$v}

# Nominal scene dimensions.
# Scene Z height had better be a multiple of ${ZFOC_STEP}.
SCENE_WX := 160
SCENE_WY := 120
SCENE_ZMIN := 30
SCENE_ZMAX := 90

ifeq "/${ZOPT}" "/S"
  # Single image at midway focus plane {Z}:
  ZFOC_MIN := ${shell echo '(${SCENE_ZMIN} + ${SCENE_ZMAX})/2.0' | bc -lq}
  ZFOC_MAX := ${ZFOC_MIN}
else ifeq "/${ZOPT}" "/M"
  # Stack of images at multiple focus plane {Z}:
  ZFOC_MIN := ${shell v=`echo '${SCENE_ZMIN} - ${ZFOC_STEP}/2.0' | bc -lq` ; printf "%.3f" $$v}
  ZFOC_MAX := ${shell v=`echo '${SCENE_ZMAX} + ${ZFOC_STEP}/2.0' | bc -lq` ; printf "%.3f" $$v}
else
  # Unknown ${ZOPT}:
  :=:
endif

PATTERN_DIR := in

SIZE_TAG := ${shell printf "%04dx%04d" ${IMAGE_NX} ${IMAGE_NY}}
SAMPLING_TAG := ${shell printf "hs%02d-kr%02d" ${PIX_SAMPLING} ${DIR_SAMPLING}}
STACK_DIR := out/st${SCENE_TYPE}-${PATTERN}-${SIZE_TAG}-${SAMPLING_TAG}

Z_DEP_FMT := ${shell printf "%08.4f" ${Z_DEP}}

SHARP_IMAGE := ${STACK_DIR}/frame-sharp/sVal.png
PIX_PLOT_FILE := ${STACK_DIR}/pixplot.txt

single-stack-clean:
	rm -fv ${STACK_DIR}/{,frame-*/}*.{ppm,pgm,png,txt}

single-stack: single-stack-clean ${SHARP_IMAGE}
	@echo "SHARP_IMAGE = ${SHARP_IMAGE}" 1>&2

${SHARP_IMAGE}: ${PROG} ${RUN_SCRIPT} Makefile
	mkdir -pv ${STACK_DIR}
	rm -f ${STACK_DIR}/selected_pixels.png
	rm -f ${STACK_DIR}/pixel-data-*.{png,txt}
	rm -f ${STACK_DIR}/frame-sharp/*.{ppm,pgm,png,txt}
	rm -f ${STACK_DIR}/frame-*-df{Z_DEP_FMT}/*.{ppm,pgm,png,txt}
	./${PROG}  \
          -imageSize  ${IMAGE_NX} ${IMAGE_NY} \
          -sceneType ${SCENE_TYPE} \
          -sceneSize 0 ${SCENE_WX}  0 ${SCENE_WY}  ${SCENE_ZMIN} ${SCENE_ZMAX} \
          -pixSampling ${PIX_SAMPLING} \
          -dirSampling ${DIR_SAMPLING} \
          -focusHeight ${ZFOC_MIN} to ${ZFOC_MAX} step ${ZFOC_STEP} \
          -dephOfFocus ${Z_DEP} \
          -patternFile ${PATTERN_DIR}/${PATTERN}.png \
          -outDir ${STACK_DIR}
          
single-stack-show-images: show_mf_images.sh ${SHARP_IMAGE}
	./show_mf_images.sh ${STACK_DIR} ${Z_DEP_FMT}
        
single-stack-plot-profiles: plot_pixel_data.sh        
	if [[ "/${ZOPT}" == "/M" ]]; then \
	  for ff in ${STACK_DIR}/pixel-data-*.txt ; do \
            plot_pixel_data.sh $${ff} ; \
          done; \
        else \
          echo "** single frame run - profiles not plotted" 1>&2; exit 1; \
        fi

endif
# End ${TARGS} section
#######################################################################
