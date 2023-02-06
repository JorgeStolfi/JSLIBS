# Last edited on 2023-02-04 07:05:52 by stolfi

PROG = minu_do_tests

TEST_LIB := libminu.a
TEST_LIB_DIR := ../..

PSVIEW := evince

JS_LIBS := \
  libeps.a \
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

