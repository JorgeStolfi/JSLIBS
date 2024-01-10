# Last edited on 2008-01-19 14:39:54 by stolfi

PROG := zftest

TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libia.a \
  libflt.a \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
