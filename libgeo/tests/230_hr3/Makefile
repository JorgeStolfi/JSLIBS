# Last edited on 2024-11-08 11:24:17 by stolfi

PROG = test_hr3
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
