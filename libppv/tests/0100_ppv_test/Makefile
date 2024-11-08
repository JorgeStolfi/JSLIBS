# Last edited on 2024-11-08 11:36:55 by stolfi

PROG = test_ppv

TEST_LIB := libppv.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
