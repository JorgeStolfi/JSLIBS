# Last edited on 2020-10-11 02:54:04 by jstolfi

PROG = test_paint

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

OUTNAME := out/img

do-test: ${PROG}
	${PROG} 
	display -filter Box -resize '400%' ${OUTNAME}-*.ppm

clean::
	/bin/rm -fv out-*.ppm
