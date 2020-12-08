# Last edited on 2008-01-14 23:07:11 by stolfi

PROG := test_order
 
TEST_LIB := libintg.a
TEST_LIB_DIR := ../..
  
JS_LIBS = \
  libgeo.a \
  libjs.a \
  libps.a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

run: ${PROG} show-error-plots
	${PROG}
	show-error-plots
 
