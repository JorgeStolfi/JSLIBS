# Last edited on 2023-01-10 16:21:17 by stolfi

PROG = waves_test

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

do-test: ${PROG}
	${PROG}
	display -title '%f' -filter Box -resize '400%' `ls out/test-*.pgm | sort`

clean::
	/bin/rm -fv out/test-*.pgm
