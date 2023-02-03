# Last edited on 2023-02-01 06:25:27 by stolfi

SHELL := /bin/bash

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := mf_0100_make_image

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

NX := 160
NY := 120

SOME_TESTS := \
  cloth-84:S \
  fiber-01:S \
  grave-02:S \
  pebbl-87:S \
  rebar-03:S \
  trees-07:S \
  wavys-09:S \
   \
  wavys-09:M \
  fiber-01:M

SOME_MORE_TESTS := \
  chalk-82:S \
  cloth-84:S \
  decor-82:S \
  falls-83:S \
  fiber-01:S \
  grave-02:S \
  india-81:S \
  pebbl-87:S \
  rebar-03:S \
  trees-07:S \
  wavys-09:S \
   \
  wavys-09:M

ALL_TESTS := \
  chalk-82:S \
  cloth-84:S \
  cloth-86:S \
  coral-03:S \
  decor-82:S \
  fallf-01:S \
  falls-83:S \
  falls-87:S \
  fiber-01:S \
  gecko-02:S \
  grass-02:S \
  grass-04:S \
  grass-09:S \
  grave-02:S \
  gravl-87:S \
  husks-01:S \
  india-81:S \
  leafy-09:S \
  leafy-82:S \
  leafy-86:S \
  mosic-96:S \
  pebbl-87:S \
  rebar-03:S \
  sqirt-01:S \
  tampi-86:S \
  tarmc-81:S \
  tarmc-82:S \
  trees-07:S \
  wavys-09:S \
  wavys-33:S \
  wheat-81:S \
   \
  wavys-09:M

check: some-tests

# SCENE_TYPES := T F R
SCENE_TYPES := F

some-tests:
	for targs in ${SOME_TESTS} ; do \
	  for st in ${SCENE_TYPES} ; do \
            ${MAKE} TARGS=$${targs} SCENE_TYPE=$${st} single-test ; \
          done ; \
        done

######################################################################
# recursive section - caller must define ${TARGS}, ${SCENE_TYPE}
ifneq "/${TARGS}" "/"
ifneq "/${SCENE_TYPE}" "/"

TEST_ARGS := ${subst :, , ${TARGS}}
PATTERN := ${word 1, ${TEST_ARGS}}
ZOPT :=    ${word 2, ${TEST_ARGS}}

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
  ZFOC_MAX := 29.5
  ZFOC_STEP := 3.0
else
  # Unknown ${ZOPT}:
  :=:
endif

PATTERN_DIR := in/pgm-256x256
PATTERN_FILE := ${PATTERN_DIR}/${PATTERN}.pgm

SIZE_TAG := ${shell printf "%04dx%04d" ${NX} ${NY}}
OUT_PREFIX := out/st${SCENE_TYPE}-${SIZE_TAG}-${PATTERN}/frame
SHARP_IMAGE := ${OUT_PREFIX}-sharp-cs.ppm

RUN_SCRIPT := run_make_image.sh

clean-out-test:
	rm -fv ${OUT_PREFIX}-*.{ppm,pgm}

single-test: clean-out-test ${SHARP_IMAGE}
	@echo "SHARP_IMAGE = ${SHARP_IMAGE}" 1>&2

${SHARP_IMAGE}: ${PROG} ${RUN_SCRIPT} Makefile
	./${RUN_SCRIPT} ${SCENE_TYPE} ${NX} ${NY} ${PATTERN_DIR} ${PATTERN} ${Z_DEP} ${ZFOC_MIN} ${ZFOC_MAX} ${ZFOC_STEP}

endif
endif
# End ${TARGS}, ${SCENE_TYPE} section
#######################################################################
