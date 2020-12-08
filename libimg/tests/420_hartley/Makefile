# Last edited on 2020-10-01 21:19:02 by jstolfi

PROG = test_hartley

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libfftw3.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

NX := 72
NY := 44
MAG := 400%

do-test: ${PROG}
	${PROG} ${NX} ${NY}
	display -title '%f' -filter Point -resize '${MAG}' `ls out-*.ppm | sort`

clean::
	/bin/rm -fv out-*.ppm
