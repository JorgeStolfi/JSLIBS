# Last edited on 2024-11-08 11:28:50 by stolfi

PROG = test_test_image

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

do-test: ${PROG}
	${PROG}
	display -title '%f' `ls out/test-*.png | sort`

clean::
	/bin/rm -fv out/test-*.ppm
