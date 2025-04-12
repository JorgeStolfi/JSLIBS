# Last edited on 2025-04-02 08:13:56 by stolfi

PROG = test_sve_step
 
TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libeigen.a \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	time ${PROG}
