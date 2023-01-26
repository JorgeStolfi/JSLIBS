# Last edited on 2023-01-24 22:10:22 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := mf_0200_op_correl

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a

OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: check clean-out generate show regress

all: check
  
DUH_PATTERNS := \
  decor-82 \
  falls-83 \
   \
  india-81 \
   \
  trees-07 \


BAD_PATTERNS := \
  chalk-82

Z_FOCUS := 15.0
FOC_DEPTH := 6.0
SCENE_TYPE := R
FD_TAG := ${shell printf "fd%08.4f" ${FOC_DEPTH}}
ZF_TAG := ${shell printf "z%08.4f" ${Z_FOCUS}}

PATTERN_SET := BEST
BASIS_NAME := DIFF
WHICH_TERMS := 1
SQR_SHARPNESS := 0
NOISE := 0.24

TERMS_TAG := pt${PATTERN_SET}-bt${BASIS_NAME}-trm${WHICH_TERMS}-qsh${SQR_SHARPNESS}-ns${NOISE}

OUT_PREFIX := out/test-st${SCENE_TYPE}-${ZF_TAG}-${TERMS_TAG}

NW := 3

TERMS_FILE := term-weights/${TERMS_TAG}.txt

TERMS_COUNT := ${shell cat ${TERMS_FILE} | egrep -e '[0-9][.]' | wc -l}

ifeq "${SQR_SHARPNESS}" "1"
  SQR_SHARPNESS_TITLE := Squared
else ifeq "${SQR_SHARPNESS}" "0"
  SQR_SHARPNESS_TITLE := Unsquared
else
  :=:
endif

ifeq "${WHICH_TERMS}" "1"
  TERMS_TITLE := All squares
else ifeq "${WHICH_TERMS}" "2"
  TERMS_TITLE := Order 2 squares
else ifeq "${WHICH_TERMS}" "3"
  TERMS_TITLE := All quadratic
else
  :=:
endif

ifeq "/${PATTERN_SET}" "/BEST"
  PATTERNS := \
    cloth-84 \
    grave-02 \
    wavys-09
else ifeq "/${PATTERN_SET}" "/SOSO"
  PATTERNS := \
    pebbl-87 \
    rebar-03 \
    fiber-01
else
  :=:
endif

IN_PREFIXES := ${addprefix in/test-st${SCENE_TYPE}-, ${addsuffix -${FD_TAG}-${ZF_TAG}, ${PATTERNS}}}

PLOT_TITLE := "Basis ${BASIS_NAME} - ${TERMS_TITLE} terms (${TERMS_COUNT}) - ${SQR_SHARPNESS_TITLE} - noise ${NOISE}"

RDATA_FILE := ${OUT_PREFIX}-rdata.txt
IMAGE_FILE := ${OUT_PREFIX}-czdsfex.ppm

check: generate show regress

clean-out:
	rm -f ${OUT_PREFIX}*.ppm ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.txt

generate: ${RDATA_FILE}

${RDATA_FILE}: ${PROG} Makefile ${TERMS_FILE}
	rm -f ${OUT_PREFIX}*.ppm ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.txt
	./${PROG} \
          ${addprefix -inPrefix , ${IN_PREFIXES}} \
          -winSize ${NW} \
          -basisName ${BASIS_NAME} \
          -termsFile ${TERMS_FILE} \
          -squared ${SQR_SHARPNESS} \
          -noise ${NOISE} \
          -outPrefix ${OUT_PREFIX}

show: ${IMAGE_FILE}

${IMAGE_FILE}: ${RDATA_FILE}
	./plot_scores_by_Z.sh ${OUT_PREFIX}
	./show_scores.sh ${word 1, ${IN_PREFIXES}} ${OUT_PREFIX}

regress: ${RDATA_FILE}
	do_regression.sh ${OUT_PREFIX}
	combine_${BASIS_NAME}_coeffs.gawk < ${OUT_PREFIX}-form.txt > .new-terms; cat .new-terms 1>&2
	@echo "should replace ${TERMS_FILE}" 1>&2
	plot_score_regression.sh ${OUT_PREFIX} ${SQR_SHARPNESS} ${PLOT_TITLE}
	plot_regression_result.sh SHOW ${OUT_PREFIX}-regr
	./compute_score_error_stats.gawk ${OUT_PREFIX}-cdata.txt


