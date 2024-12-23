# Last edited on 2023-10-05 06:30:13 by stolfi

PROG := test_haf
 
TEST_LIB := libquad.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}

