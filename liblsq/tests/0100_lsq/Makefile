# Last edited on 2024-12-21 05:06:32 by stolfi

PROG = test_lsq
 
TEST_LIB := liblsq.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	time ${PROG}
