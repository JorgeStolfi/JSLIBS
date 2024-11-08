# Last edited on 2024-11-08 12:30:24 by stolfi

PROG := test_fp

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

run: ${PROG}
	${PROG}
