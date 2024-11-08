# Last edited on 2024-11-04 06:56:53 by stolfi

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
	rm -f out/*.eps
	${PROG}
	for f in out/*opt01*.eps ; do \
          echo $$f; ${PSVIEW} $${f/opt01/opt*} ; \
        done

clean::
	-/bin/rm -f *.ps *.eps

