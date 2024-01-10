# Last edited on 2014-07-27 03:05:19 by stolfilocal

PROG = test_hermite3

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: clean do-test

NAME := ta

do-test: ${PROG}
	mkdir -p out
	${PROG} ${NAME}
	for f in `ls out/${NAME}*.txt | sort` ; do \
          show-plot.sh $$f ; \
        done

clean::
	/bin/rm -fv out/*.txt out/*.eps
