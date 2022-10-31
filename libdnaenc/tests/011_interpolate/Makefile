# Last edited on 2022-10-31 14:37:34 by stolfi

TEST_LIB := libdnaenc.a
TEST_LIB_DIR := ../..
PROG := test_dnae_interpolate

JS_LIBS := \
  libmsmatch.a \
  libgeo.a \
  libeps.a \
  libimg.a \
  libjs.a

OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

STEPS := 0.1250 0.2500 0.5000 1.0000

all: check

check:  clean run-all

run-all: 
	for step in ${STEPS} ; do \
            ${MAKE} STEP=$$step clean-single run-single show-single ; \
        done
     
# ======================================================================
# Recursive make - must define ${STEP}
ifneq "/${STEP}" "/"        

IN_DIR := ${STOLFIHOME}/programs/c/DNA/DATA/small/eqs

SEQ := 100A-wg020-wg060-02
XSTEP := ${subst .,,${STEP}}

IN_NAME := ${IN_DIR}/${SEQ}

OUT_NAME := out/${SEQ}-${XSTEP}

SHOWPS := evince
  
run-single: ${PROG}
	${PROG} \
          -step ${STEP} \
          ${IN_NAME} \
          ${OUT_NAME}

EPS_FILE := ${OUT_NAME}-ot.eps

show-single: ${EPS_FILE}
	${SHOWPS} ${EPS_FILE}
  
clean-single:: 
	-/bin/rm -f ${OUT_NAME}-*.pgm ;
	-/bin/rm -f ${OUT_NAME}-*.egs ;
	-/bin/rm -f ${OUT_NAME}-*.eps ;
 
endif
# ======================================================================
