# Last edited on 2024-12-05 22:26:18 by stolfi

PROG = test_paint

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libgausol.a \
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
