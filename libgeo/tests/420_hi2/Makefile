# Last edited on 2024-11-08 12:46:31 by stolfi

PROG = test_hi2
 
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
