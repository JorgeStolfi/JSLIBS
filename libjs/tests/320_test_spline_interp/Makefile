# Last edited on 2013-10-28 20:29:38 by stolfilocal

PROG = test_spline_interp

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
	${PROG} ${NAME}
	for f in `ls out/${NAME}-r*-o*-k*.txt | sort` ; do \
          show-plot.sh $$f ; \
        done

clean::
	/bin/rm -fv out/*.txt out/*.eps
