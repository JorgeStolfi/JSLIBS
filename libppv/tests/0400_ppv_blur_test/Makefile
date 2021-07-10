# Last edited on 2021-07-03 21:33:30 by jstolfi

PROG = ppv_blur_test

TEST_LIB := libppv.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libimg.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so  

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  do-test

do-test: ${PROG}
	rm -fv out/*.png
	${PROG}
	eom out/*.png
