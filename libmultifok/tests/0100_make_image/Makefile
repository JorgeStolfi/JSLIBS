# Last edited on 2024-10-01 16:10:42 by stolfi

SHELL := /bin/bash

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_make_image

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: \
  some-tests \
  single-test \
  clean_out-test

all: check

# Image dimensions:
# IMAGE_NX := 640
# IMAGE_NY := 480

# IMAGE_NX := 160
# IMAGE_NY := 120

IMAGE_NX := 320
IMAGE_NY := 240

# Nominal scene dimensions:
SCENE_WX := 160
SCENE_WY := 120

# {TEXTURE}:{MULTIPLICITY}:{SCENE_TYPE}
# where {MULTIPLICITY} is "S" for single frame, "M" for frame stack,
# and {SCENE_TYPE} is "R" for ramp only, "F" for non-overlapping objs,
# "T" for overlapping.

QUICK_TESTS := \
  wavys-09:M:F \
  cloth-84:S:F

SOME_TESTS := \
  cloth-84:S:R \
  fiber-01:S:R \
  fiber-01:S:R \
  grave-02:S:R \
  pebbl-87:S:R \
  wavys-09:S:R \
  wavyu-zz:S:R \
   \
  cloth-84:S:F \
  grave-02:S:F \
  wavys-09:S:F \
   \
  wavys-09:M:F \
  fiber-01:M:F \

ALL_TESTS := \
  cloth-84:S:R \
  fiber-01:S:R \
  fiber-01:S:R \
  grave-02:S:R \
  pebbl-87:S:R \
  wavys-09:S:R \
  wavyu-zz:S:R \
   \
  chalk-82:S:F \
  cloth-84:S:F \
  cloth-86:S:F \
  coral-03:S:F \
  decor-82:S:F \
  fallf-01:S:F \
  falls-83:S:F \
  falls-87:S:F \
  fiber-01:S:F \
  gecko-02:S:F \
  grass-02:S:F \
  grass-04:S:F \
  grass-09:S:F \
  grave-02:S:F \
  gravl-87:S:F \
  husks-01:S:F \
  india-81:S:F \
  leafy-09:S:F \
  leafy-82:S:F \
  leafy-86:S:F \
  mosic-96:S:F \
  pebbl-87:S:F \
  rebar-03:S:F \
  sqirt-01:S:F \
  tampi-86:S:F \
  tarmc-81:S:F \
  tarmc-82:S:F \
  trees-07:S:F \
  wavys-09:S:F \
  wavys-33:S:F \
  wheat-81:S:F \
   \
  wavys-09:M:F \
  fiber-01:M:F

TESTS := ${QUICK_TESTS}

check: 
	for targs in ${TESTS} ; do \
	  ${MAKE} TARGS=$${targs} single-test ; \
        done

######################################################################
# recursive section - caller must define ${TARGS}
ifneq "/${TARGS}" "/"

TEST_ARGS :=  ${subst :, , ${TARGS}}
PATTERN :=    ${word 1, ${TEST_ARGS}}
ZOPT :=       ${word 2, ${TEST_ARGS}}
SCENE_TYPE := ${word 3, ${TEST_ARGS}}

ifeq "/${ZOPT}" "/S"
  # Single image at midway focus plane {Z}:
  Z_DEP := 6.0
  ZFOC_MIN := 15.5
  ZFOC_MAX := 15.5
  ZFOC_STEP := 3.0
else ifeq "/${ZOPT}" "/M"
  # Stack of images at multiple focus plane {Z}:
  Z_DEP := 6.0
  ZFOC_MIN := 0.5
  ZFOC_MAX := 30.5
  ZFOC_STEP := 3.0
else
  # Unknown ${ZOPT}:
  :=:
endif

PATTERN_DIR := in/pgm-256x256
PATTERN_FILE := ${PATTERN_DIR}/${PATTERN}.pgm

SIZE_TAG := ${shell printf "%04dx%04d" ${IMAGE_NX} ${IMAGE_NY}}
OUT_PREFIX := out/st${SCENE_TYPE}-${SIZE_TAG}-${PATTERN}/frame
SHARP_IMAGE := ${OUT_PREFIX}-sharp-cs.ppm
PIX_PLOT_FILE := ${OUT_PREFIX}-pixplot.txt

RUN_SCRIPT := run_mfok_make_image.sh

clean-out-test:
	rm -fv ${OUT_PREFIX}-*.{ppm,pgm}

single-test: clean-out-test ${SHARP_IMAGE}
	@echo "SHARP_IMAGE = ${SHARP_IMAGE}" 1>&2

${SHARP_IMAGE}: ${PROG} ${RUN_SCRIPT} Makefile
	./${RUN_SCRIPT} \
          ${SCENE_TYPE} ${SCENE_WX} ${SCENE_WY} \
          ${IMAGE_NX} ${IMAGE_NY} \
          ${PATTERN_DIR} ${PATTERN} \
          ${Z_DEP} ${ZFOC_MIN} ${ZFOC_MAX} ${ZFOC_STEP}
	if [[ "/${ZOPT}" == "/M" ]]; then \
          plot_hrad2_by_z.sh ${PIX_PLOT_FILE} ; \
        fi

endif
# End ${TARGS} section
#######################################################################
