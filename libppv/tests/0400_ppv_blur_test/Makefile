# Last edited on 2024-11-08 11:37:55 by stolfi

PROG = test_ppv_blur

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
