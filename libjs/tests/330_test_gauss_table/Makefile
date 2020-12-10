# Last edited on 2020-12-09 23:54:28 by jstolfi

PROG = test_gauss_table

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
