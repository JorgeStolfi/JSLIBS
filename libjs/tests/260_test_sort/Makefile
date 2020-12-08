# Last edited on 2008-01-14 23:26:28 by stolfi

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
