# Last edited on 2011-09-20 01:19:38 by stolfilocal

PROG = test_float_array

TEST_LIB := libjsarray.a
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
	rm -f out/*.txt
