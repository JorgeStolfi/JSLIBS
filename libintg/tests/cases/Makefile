# Last edited on 2023-02-03 23:24:22 by stolfi

PROG = test_cases
 
TEST_LIB := libintg.a
TEST_LIB_DIR := ../..
  
JS_LIBS = \
  libeps.a \
  libgeo.a \
  libjs.a
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

clean::
	-/bin/rm -f *.eps *.ps

run: ${PROG}
	${PROG}
	show-intg-plots
	show-steps euler 30
	show-steps rkfo4 30
  
