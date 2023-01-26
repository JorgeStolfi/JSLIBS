# Last edited on 2023-01-25 23:46:52 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := mf_0180_diff_ops

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
  india-81 \
  trees-07 \

BAD_PATTERNS := \
  chalk-82

Z_FOCUS := 15.5
FOC_DEPTH := 6.0
SCENE_TYPE := R
FD_TAG := ${shell printf "fd%05.2f" ${FOC_DEPTH}}
ZF_TAG := ${shell printf "z%04.1f" ${Z_FOCUS}}

PATTERN_SET := SOSO
BASIS_NAME := LAPL
NOISE := 0.04
UNIT_TERM := 0

PAT_BASIS_TAG := pt${PATTERN_SET}-bt${BASIS_NAME}-ns${NOISE}

IN_PREFIX := in/single-${ZF_TAG}/

OUT_PREFIX := out/test-st${SCENE_TYPE}-${ZF_TAG}-${PAT_BASIS_TAG}

NW := 3

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

IN_IMAGES := ${addprefix st${SCENE_TYPE}-, ${addsuffix -${FD_TAG}-${ZF_TAG}, ${PATTERNS}}}

ifneq "/${UNIT_TERM}" "/0"
  UNIT_TERM_TITLE := +UNIT
else
  UNIT_TERM_TITLE :=
endif

COEFF_PLOT_TITLE := "Patset ${PATTERN_SET} - ${BASIS_NAME} - Noise ${NOISE}"
REGR_PLOT_TITLE := "Patset ${PATTERN_SET} - ${BASIS_NAME}${UNIT_TERM_TITLE} - Noise ${NOISE}"

ODATA_FILE := ${OUT_PREFIX}-odata.txt
IMAGE_FILE := ${OUT_PREFIX}-czdsfex.ppm

check: generate show regress

clean-out:
	rm -f ${OUT_PREFIX}*.ppm ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.txt

generate: ${ODATA_FILE}

${ODATA_FILE}: ${PROG} Makefile
	rm -f ${OUT_PREFIX}*.ppm ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.txt
	./${PROG} \
          -inPrefix ${IN_PREFIX} \
          ${addprefix -image , ${IN_IMAGES}} \
          -basisName ${BASIS_NAME} \
          -noise ${NOISE} \
          -outPrefix ${OUT_PREFIX}

show: ${IMAGE_FILE}

${IMAGE_FILE}: ${ODATA_FILE}
	./plot_ops_squared_by_sharp.sh ${OUT_PREFIX} ${BASIS_NAME} ${UNIT_TERM} ${COEFF_PLOT_TITLE}
	./show_ops_squared.sh ${IN_PREFIX}${word 1,${IN_IMAGES}} ${OUT_PREFIX}
 
# regress: ${ODATA_FILE}
# 	do_regression.sh ${OUT_PREFIX}
# 	combine_${BASIS_NAME}_coeffs.gawk < ${OUT_PREFIX}-form.txt > .new-terms; cat .new-terms 1>&2
# 	@echo "should replace ${TERMS_FILE}" 1>&2
# 	plot_score_regression.sh ${OUT_PREFIX} ${SQR_SHARPNESS} ${REGR_PLOT_TITLE}
# 	plot_regression_result.sh SHOW ${OUT_PREFIX}-regr
# 	./compute_score_error_stats.gawk ${OUT_PREFIX}-cdata.txt


