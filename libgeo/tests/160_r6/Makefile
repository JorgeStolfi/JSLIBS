# Last edited on 2016-04-03 13:20:04 by stolfilocal

PROG = r6test
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

build-test: ${PROG}

do-test: ${PROG}
	${PROG}
