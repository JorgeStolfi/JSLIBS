# Last edited on 2023-02-03 23:24:45 by stolfi

PROG := test_bz

TEST_LIB := libbezier.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libeps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

clean::
	-/bin/rm -f *.eps *.ps

run: ${PROG}
	time ${PROG}
