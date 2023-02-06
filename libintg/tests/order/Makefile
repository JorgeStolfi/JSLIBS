# Last edited on 2023-02-03 23:24:15 by stolfi

PROG := test_order
 
TEST_LIB := libintg.a
TEST_LIB_DIR := ../..
  
JS_LIBS = \
  libeps.a \
  libgeo.a \
  libjs.a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  clean run

run: ${PROG} show-error-plots
	${PROG}
	show-error-plots
 
