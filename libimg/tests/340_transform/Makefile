# Last edited on 2024-12-20 18:27:06 by stolfi

PROG = test_transform

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a  \
  libgausol.a\
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
