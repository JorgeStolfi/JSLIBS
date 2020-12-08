# Last edited on 2016-12-28 23:59:15 by stolfilocal

PROG := iavalid

TEST_LIB := libia.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libflt.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

run: ${PROG}
	${PROG}
