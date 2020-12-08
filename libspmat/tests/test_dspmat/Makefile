# Last edited on 2009-01-17 19:20:28 by stolfi

PROG = test_dspmat

TEST_LIB := libspmat.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
