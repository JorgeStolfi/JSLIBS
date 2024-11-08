# Last edited on 2024-11-08 12:46:56 by stolfi

PROG := test_zf

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
