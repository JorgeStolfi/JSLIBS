# Last edited on 2024-06-22 20:54:56 by stolfi

PROG := zftest

TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libia.a \
  libflt.a \
  libeps.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
