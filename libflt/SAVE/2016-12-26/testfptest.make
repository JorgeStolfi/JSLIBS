# Last edited on 2016-12-26 22:00:51 by stolfilocal

PROG := testfptest

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

OTHER_I_FLAGS := \
  -IJSLIBS/libjs

OTHER_LIBS := \
  JSLIBS/libjs/libjs.a

include GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
