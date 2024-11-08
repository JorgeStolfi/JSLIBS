# Last edited on 2024-11-08 11:26:17 by stolfi

PROG = test_ellipse
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	time ${PROG}
