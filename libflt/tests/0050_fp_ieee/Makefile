# Last edited on 2024-11-08 12:39:36 by stolfi

PROG := test_fp_ieee

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

run: ${PROG}
	${PROG}
