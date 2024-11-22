# Last edited on 2024-11-17 12:08:43 by stolfi

PROG = test_sort
 
TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

NELEMS := 1000

do-test: ${PROG}
	time ${PROG} ${NELEMS}
