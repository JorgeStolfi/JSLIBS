# Last edited on 2009-03-12 00:06:54 by stolfi

PROG := test_sort_speed

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: do-test

do-test: ${PROG}
	time ${PROG}

OUT_NAME := ${PROG}

CMPS_FILE := ${OUT_NAME}-ncmp.tex
TIME_FILE := ${OUT_NAME}-time.tex

clean::
	-/bin/rm -f ${OUT_NAME}-ncmp.tex ${OUT_NAME}-time.tex
