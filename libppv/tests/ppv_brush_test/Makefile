# Last edited on 2021-06-13 13:46:06 by jstolfi

PROG = ppv_brush_test

TEST_LIB := libppv.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a\
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
