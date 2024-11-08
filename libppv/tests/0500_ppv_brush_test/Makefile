# Last edited on 2024-11-08 11:38:12 by stolfi

PROG = test_ppv_brush

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
