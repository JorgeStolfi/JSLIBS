# Last edited on 2024-06-28 02:11:58 by stolfi

PROG = test_jsdebug

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test show-plots

all: ${PROG}
	mkdir -p out
	${PROG}

clean::
