# Last edited on 2024-11-29 22:21:05 by stolfi

PROG = test_qmin_simplex
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a \
  libgausol.a \
  libgeo.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	time ${PROG}

clean::
