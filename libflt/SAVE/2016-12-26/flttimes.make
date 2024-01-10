# Last edited on 2016-12-26 21:44:56 by stolfilocal

PROG := flttimes

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

OTHER_I_FLAGS := \
  -IJSLIBS/libjs \
  -IJSLIBS/libflt

OTHER_LIBS := \
  JSLIBS/libjs/libjs.a

include GENERIC-LIB-TEST.make

.PHONY:: run

all: clean run

run: ${PROG}
	${PROG}
