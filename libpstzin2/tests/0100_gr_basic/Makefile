# Last edited on 2025-03-14 05:55:23 by stolfi

PROG = test_gr_basic

TEST_LIB := libpstzin2.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
  libeps.a \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a \
  libjsextra.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

do-test: ${PROG}
	mkdir -p out
	${PROG}
 
clean::
	/bin/rm -fv out/*-*.fni
