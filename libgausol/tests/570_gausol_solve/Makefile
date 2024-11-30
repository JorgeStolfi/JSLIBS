# Last edited on 2024-11-28 21:03:56 by stolfi

PROG = test_gausol_solve
 
TEST_LIB := libgausol.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}

clean::
