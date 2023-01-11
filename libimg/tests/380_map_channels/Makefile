# Last edited on 2023-01-07 13:54:32 by stolfi

PROG = test_map_channels

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

do-test: ${PROG}
	mkdir -p out
	rm -f out/*.fni out/*.ppm out/*.pgm
	${PROG}
	display -title '%f' -filter Point -resize '400%' `ls out/test-*.p?m | sort`

clean::
	/bin/rm -fv out/test-*.ppm
