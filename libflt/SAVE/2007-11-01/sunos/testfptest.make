# Last edited on 2008-01-14 22:47:01 by stolfi

PROG := testfptest

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
