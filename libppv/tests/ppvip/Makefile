# Last edited on 2008-01-14 22:47:01 by stolfi

PROG := ppvip

TEST_LIB := libppv.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  do-test

do-test: ${PROG}
	${PROG}

