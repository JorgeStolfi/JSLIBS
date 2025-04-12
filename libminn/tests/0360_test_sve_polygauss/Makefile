# Last edited on 2025-03-19 12:51:47 by stolfi

PROG = test_sve_polygauss
 
TEST_LIB := libminn.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a
  
OTHER_LIBS :=
 
TESTS :=
 
.PHONY:: all do-test clean

all: do-test
	
do-test: clean ${PROG} 
	${PROG}
	plot_tests.sh

clean::
	-rm -f out/train_*.txt out/train_*.png 

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
