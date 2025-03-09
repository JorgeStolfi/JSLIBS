# Last edited on 2025-02-19 19:16:23 by stolfi

PROG := test_interpolate

TEST_LIB := libpst.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
  libgeo.a \
  libgausol.a \
  libimg.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: run_tests.sh ${PROG}
	run_tests.sh
