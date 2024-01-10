# Last edited on 2013-12-17 00:44:07 by stolfilocal

PROG = test_fget

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: do-test

do-test: ${PROG} out
	${PROG}

out:
	mkdir -p out
