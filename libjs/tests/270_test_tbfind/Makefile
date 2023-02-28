# Last edited on 2023-02-27 06:54:52 by stolfi

PROG := test_tbfind
 
TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
