# Last edited on 2015-12-21 16:54:58 by stolfilocal

PROG := g3map_test_00_basic
 
TEST_LIB := libgmap.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libgem.a \
  libjs.a
 
all: check

check:  do-test

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

do-test: ${PROG}
	${PROG}
