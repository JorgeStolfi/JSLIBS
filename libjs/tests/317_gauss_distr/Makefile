# Last edited on 2019-12-04 19:49:25 by jstolfi

PROG = test_gauss_distr

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

NAME := ta

do-test: ${PROG}
	mkdir -p out
	${PROG}

clean::
	/bin/rm -fv out/*.txt out/*.eps
