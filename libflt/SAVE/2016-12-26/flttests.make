# Last edited on 2016-12-26 21:43:08 by stolfilocal

PROG := flttests

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

OTHER_I_FLAGS := \
  -IJSLIBS/libjs \
  -IJSLIBS/libflt

OTHER_LIBS := \
  JSLIBS/libjs/libjs.a

include GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
