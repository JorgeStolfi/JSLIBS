# Last edited on 2011-12-22 15:41:26 by stolfilocal

PROG := test_quad
 
TEST_LIB := libquad.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}

