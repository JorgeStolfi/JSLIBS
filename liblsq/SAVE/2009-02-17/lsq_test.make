# Last edited on 2008-01-14 23:22:54 by stolfi

PROG = lsq_test
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: do-test

do-test: ${PROG}
	time ${PROG}
