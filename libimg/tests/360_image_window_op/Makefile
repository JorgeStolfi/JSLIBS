# Last edited on 2020-11-15 14:04:38 by jstolfi

PROG = test_image_window_op

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
	display -title '%f' -filter Point -resize '400%' `ls out/t-*.ppm | sort`

clean::
	/bin/rm -fv out/t-*.ppm
