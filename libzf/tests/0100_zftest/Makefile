# Last edited on 2023-02-18 03:21:13 by stolfi

PROG := zftest

TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libia.a \
  libflt.a \
  libeps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
