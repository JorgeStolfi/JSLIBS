# Last edited on 2008-01-14 23:25:55 by stolfi

PROG = test_jsmath

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: do-test

do-test: ${PROG}
	${PROG}
