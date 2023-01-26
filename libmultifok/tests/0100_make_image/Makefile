# Last edited on 2023-01-25 22:55:08 by stolfi

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
  generate show export \
  single-test \
  clean_out-test generate-test show-test export_test

all: check

SCENE_TYPE := F

SOME_TESTS := \
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

ifeq "${SCENE_TYPE}" "F"
  # Non-overlapping disks/balls:
  BACK_TILT := F
  OVERLAP := F
else ifeq "${SCENE_TYPE}" "T"
  # Overlapping disks/balls:
  BACK_TILT := F
  OVERLAP := T
else ifeq "${SCENE_TYPE}" "R"
  # Tilted backplane only:
  BACK_TILT := T
  OVERLAP := F
else
  := ERROR
endif

check: generate show export

generate:
	for targs in ${SOME_TESTS} ; do \
          ${MAKE} TARGS=$${targs} generate-test ; \
        done

show:
	for targs in ${SOME_TESTS} ; do \
          ${MAKE} TARGS=$${targs} show-test ; \
        done

export:
	for targs in ${SOME_TESTS} ; do \
  	  ${MAKE} TARGS=$${targs} export-test ; \
        done

######################################################################
# recursive section - caller must define ${TARGS}
ifneq "/${TARGS}" "/"

TEST_ARGS := ${subst :, , ${TARGS}}
PATTERN := ${word 1, ${TEST_ARGS}}
ZOPT :=    ${word 2, ${TEST_ARGS}}

ifeq "/${ZOPT}" "/S"
  # Single image at midway focus plane {Z}:
  ZFOC_MIN := 15.5
  ZFOC_MAX := 15.5
  ZFOC_STEP := 3.0
  ZFOCS := -zRange 15.5 15.5 -zStep 1.0
  OUT_SUBDIR := single-z${shell printf "%04.1f" "${ZFOC_MIN}"}
else ifeq "/${ZOPT}" "/M"
  # Stack of images at multiple focus plane {Z}:
  ZFOC_MIN := 0.5
  ZFOC_MAX := 29.5
  ZFOC_STEP := 3.0
  OUT_SUBDIR := stack-z${shell printf "%04.1f-%04.1f-%04.1f" "${ZFOC_MIN}" "${ZFOC_MAX}" "${ZFOC_STEP}" }
else
  # Unknown ${ZOPT}:
  :=:
endif

NX := 160
NY := 120
NP := 5
NR := 40
FOC_DEPTH := 6.0

PATTERN_FILE := in/pgm-256x256/${PATTERN}.pgm

OUT_PREFIX := out/${OUT_SUBDIR}/st${SCENE_TYPE}-${PATTERN}

SHARP_IMAGE := ${OUT_PREFIX}-sharp-c.ppm

clean_out-test:
	mkdir -p out/${OUT_SUBDIR}
	rm -f ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.ppm

generate-test: ${SHARP_IMAGE} clean_out-test

${SHARP_IMAGE}: ${PROG}
	mkdir -p out/${OUT_SUBDIR}
	./${PROG} \
          -imgSize ${NX} ${NY} \
          -backTilt ${BACK_TILT} \
          -overlap ${OVERLAP} \
          -pixSampling ${NP} \
          -dirSampling ${NR} \
          -zRange ${ZFOC_MIN} ${ZFOC_MAX} -zStep ${ZFOC_STEP} \
          -focDepth ${FOC_DEPTH} \
          -pattern ${PATTERN_FILE} \
          -outPrefix ${OUT_PREFIX}

show-test: ${SHARP_IMAGE}
	./show_images.sh ${OUT_PREFIX}

export-test: ${SHARP_IMAGE}
	@echo "!! Not exported" 1>&2

endif
# End ${TARGS} section
#######################################################################
