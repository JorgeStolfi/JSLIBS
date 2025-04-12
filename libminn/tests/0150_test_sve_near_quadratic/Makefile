# Last edited on 2025-03-19 12:51:12 by stolfi

PROG = test_sve_near_quadratic
 
TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libgausol.a \
  libjs.a
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: 

all: check

check:  ${PROG} 
	time ${PROG}

clean::
	-/bin/rm -f *-ot.dat
