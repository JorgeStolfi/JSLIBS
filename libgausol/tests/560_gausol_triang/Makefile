# Last edited on 2024-11-28 15:31:06 by stolfi

PROG = test_gausol_triang 
 
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
