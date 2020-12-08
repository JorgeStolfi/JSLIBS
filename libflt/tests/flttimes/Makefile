# Last edited on 2008-01-15 02:55:25 by stolfi

PROG := flttimes

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
