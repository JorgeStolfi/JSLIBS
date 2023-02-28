# Last edited on 2023-02-26 08:53:15 by stolfi

PROG := testfig

TEST_LIB := libeps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test show-epsf

all: check

check:  do-test

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/test_fig

MAIN_EPS := ${OUT_PREFIX}_00000_A.eps

# PSVIEW := ghostview
# PSVIEW := gv
PSVIEW := evince
        
do-test: ${MAIN_EPS} show-epsf
  
${MAIN_EPS} : ${PROG} ${OUT_DIR}
	${PROG}

show-epsf: ${MAIN_EPS}
	for f in `ls ${OUT_PREFIX}_*.eps | sort` ; do \
	  ${PSVIEW} $$f ; \
        done
        
${OUT_DIR}:
	mkdir -p ${OUT_DIR}

clean::
	-/bin/rm -f ${OUT_PREFIX}_*.eps
  
