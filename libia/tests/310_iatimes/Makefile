# Last edited on 2016-12-28 23:58:51 by stolfilocal

PROG := iatimes

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
