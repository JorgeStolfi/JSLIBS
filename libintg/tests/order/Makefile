# Last edited on 2023-03-29 19:37:10 by stolfi

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

run: ${PROG} show_error_plots.sh
	${PROG}
	show_error_plots.sh
 
