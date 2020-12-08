# Last edited on 2017-04-17 20:19:21 by stolfilocal

PROG = test_image_test

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
	display -title '%f' `ls out/test-*.ppm | sort`

clean::
	/bin/rm -fv out/test-*.ppm
