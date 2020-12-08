# Last edited on 2008-01-14 22:43:17 by stolfi

PROG := aatests

TEST_LIB := libaa.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libia.a \
  libflt.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

run: ${PROG}
	${PROG}
