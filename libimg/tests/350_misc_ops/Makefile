# Last edited on 2020-10-11 15:27:22 by jstolfi

PROG = test_misc_ops

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

RUNS := \
  thin:pgm:D \
  
MORE_RUNS :=  \
  lite:ppm:G \
  lite:ppm:R \
  dark:ppm:G \
  dark:ppm:R

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-all run-single show-output plot-output

all: check

check:  clean run-all

run-all:  
	for r in ${RUNS}; do \
          ${MAKE} RUN=$$r run-single; \
        done

OUT_DIR := out

OUTNAME := ${OUT_DIR}/img

        
clean::
	/bin/rm -fv ${OUTNAME}-*.pgm ${OUTNAME}-*.ppm

########################################################################        
# Section for recursive "make" -- client must defined ${RUN}

ifneq "/${RUN}" "/"

RUN_FIELDS := ${subst :, ,${RUN}}
NAME :=    ${word 1, ${RUN_FIELDS}}
EXT :=     ${word 2, ${RUN_FIELDS}}
OP_TAG :=  ${word 3, ${RUN_FIELDS}}

ifeq "/${OP_TAG}" "/G"
  OP_OPTION := -gradSqr
endif
ifeq "/${OP_TAG}" "/R"
  OP_OPTION := -gradSqrRel
endif
ifeq "/${OP_TAG}" "/D"
  OP_OPTION := -dilate
endif

IN_FILE := in/${NAME}.${EXT}
OUT_FILE := ${OUT_DIR}/${NAME}-${OP_TAG}.${EXT}

run-single: ${OUT_FILE} show-output plot-output

${OUT_FILE}: ${PROG} ${OUT_DIR} ${IN_FILE}
	${PROG} ${OP_OPTION} < ${IN_FILE} > ${OUT_FILE}

show-output:
	display -title '%f' -gamma 1.0 ${OUT_FILE} ${IN_FILE}
        
TEMP := /tmp/tgr

${OUT_DIR}:
	mkdir -p ${OUT_DIR}

plot-output:
	pnm_to_fni \
            -min 0 -max 1 \
          < ${OUT_FILE} \
          > ${TEMP}.fni
	fni_plot.sh \
            -title "${OUT_FILE}" \
            -range -1.01 2.01 \
          < ${TEMP}.fni \
          > ${TEMP}.eps
	-display -title '%f' -background White +matte ${TEMP}.eps

endif
# End ${RUN} section
########################################################################
