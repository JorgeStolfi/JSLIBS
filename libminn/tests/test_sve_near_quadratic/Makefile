# Last edited on 2009-02-22 19:20:11 by stolfi

PROG = test_sve_near_quadratic
 
TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: 

all: check

check:  ${PROG} 
	time ${PROG}

clean::
	-/bin/rm -f *-ot.dat
