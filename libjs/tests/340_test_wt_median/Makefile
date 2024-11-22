# Last edited on 2024-11-19 06:46:24 by stolfi

PROG = test_wt_median

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
