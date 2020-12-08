# Last edited on 2011-05-29 10:51:20 by stolfi

PROG = minu_do_tests

TEST_LIB := libminu.a
TEST_LIB_DIR := ../..

PSVIEW := okular

JS_LIBS := \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	rm -f *.ps *.eps
	${PROG}
	for f in *.ps *.eps ; do \
          echo $$f; ${PSVIEW} $$f ; \
        done

clean::
	-/bin/rm -f *.ps *.eps

