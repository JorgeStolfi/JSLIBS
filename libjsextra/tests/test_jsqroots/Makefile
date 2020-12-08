# Last edited on 2008-07-24 18:52:10 by stolfi

PROG = testjsqroots

TEST_LIB := libjsextra.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}
