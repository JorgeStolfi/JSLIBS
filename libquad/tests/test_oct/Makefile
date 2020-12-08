# Last edited on 2009-03-07 16:34:12 by stolfi

PROG := test_oct
 
TEST_LIB := libquad.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}

