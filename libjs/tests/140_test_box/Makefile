# Last edited on 2021-12-31 21:01:38 by stolfi

PROG := test_box
 
TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
