# Last edited on 2024-11-08 12:34:18 by stolfi

PROG := test_flt

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

run: ${PROG}
	${PROG}
