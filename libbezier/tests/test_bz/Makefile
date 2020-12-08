# Last edited on 2008-01-14 22:45:47 by stolfi

PROG := test_bz

TEST_LIB := libbezier.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

clean::
	-/bin/rm -f *.eps *.ps

run: ${PROG}
	time ${PROG}
