# Last edited on 2016-12-26 22:32:18 by stolfilocal

PROG := zftest

TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

OTHER_I_FLAGS := \
  -IJSLIBS/libjs \
  -IJSLIBS/libps \
  -IJSLIBS/libflt \
  -IJSLIBS/libia \
  -IJSLIBS/libaa \
  -IJSLIBS/libfgraph

OTHER_LIBS := \
  JSLIBS/libfgraph/libfgraph.a \
  JSLIBS/libaa/libaa.a \
  JSLIBS/libia/libia.a \
  JSLIBS/libflt/libflt.a \
  JSLIBS/libps/libps.a \
  JSLIBS/libjs/libjs.a

include GENERIC-LIB-TEST.make

all: clean run

run: ${PROG}
	${PROG}
