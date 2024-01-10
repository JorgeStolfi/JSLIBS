# Last edited on 2012-02-25 07:38:11 by stolfilocal

PROG = test_jsrandom

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: do-test

do-test: ${PROG}
	${PROG}
	plot-hist.sh uint64_random-lo-A
	plot-hist.sh uint64_random-hi-A
	plot-hist.sh abrandom-A
	plot-hist.sh abrandom-B
	plot-hist.sh drandom-A

clean::
	rm -f out/*.his out/*.png
