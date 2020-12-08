# Last edited on 2011-06-06 22:30:39 by stolfi

PROG = test_now
 
TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

OUT_NAME := out/${PROG}

do-test: ${PROG}
	@echo "gnuplot did not run" > ${OUT_NAME}.err
	${PROG}
	@cat ${OUT_NAME}.err

clean:: 
	-/bin/rm -f ${OUT_NAME}.dat ${OUT_NAME}.gpl ${OUT_NAME}.err 

