# Last edited on 2025-01-19 23:36:40 by stolfi

PROG := make_test_maps

TEST_LIB := libpst.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
  libgeo.a \
  libgausol.a \
  libimg.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: run_tests.sh ${PROG}
	run_tests.sh
        

