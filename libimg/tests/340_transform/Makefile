# Last edited on 2012-02-20 05:58:41 by stolfilocal

PROG = test_transform

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
	display -title '%f' -filter Point -resize '400%' `ls out/test-*.ppm | sort`

clean::
	/bin/rm -fv out/test-*.ppm
