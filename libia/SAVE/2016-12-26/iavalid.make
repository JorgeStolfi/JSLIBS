# Last edited on 2016-12-26 22:03:30 by stolfilocal

PROG := iavalid

TEST_LIB := libia.a
TEST_LIB_DIR := ../..

OTHER_I_FLAGS := \
  -IJSLIBS/libjs \
  -IJSLIBS/libflt

OTHER_LIBS := \
  JSLIBS/libflt/libflt.a \
  JSLIBS/libjs/libjs.a

include GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
