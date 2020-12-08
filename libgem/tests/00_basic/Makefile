# Last edited on 2014-07-11 00:24:46 by stolfilocal

PROG := gem_test_00_basic
 
TEST_LIB := libgem.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}
