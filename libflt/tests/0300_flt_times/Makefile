# Last edited on 2024-11-08 12:36:55 by stolfi

PROG := test_flt_times

TEST_LIB := libflt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run

all: check

check:  clean run

run: ${PROG}
	${PROG}
