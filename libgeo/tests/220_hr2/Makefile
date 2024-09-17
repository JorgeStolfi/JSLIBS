# Last edited on 2024-09-17 17:29:59 by stolfi

PROG = hr2_test
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
