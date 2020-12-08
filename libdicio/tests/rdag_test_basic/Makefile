# Last edited on 2009-10-25 16:26:53 by stolfi

PROG = rdag_test_basic
 
TEST_LIB := libdicio.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
