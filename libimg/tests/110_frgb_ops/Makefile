# Last edited on 2010-05-24 23:26:45 by stolfi

PROG = test_frgb_ops

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-tests 

all: check

check:  clean run-tests

run-tests: ${PROG}
	${PROG}

        
clean::
