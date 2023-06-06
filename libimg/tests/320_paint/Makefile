# Last edited on 2023-04-23 11:29:39 by stolfi

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
	display -filter Box -resize '200%' ${OUTNAME}-*.ppm
	display ${OUTNAME}-*.ppm

clean::
	/bin/rm -fv out-*.ppm
