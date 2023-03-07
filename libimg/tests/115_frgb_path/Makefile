# Last edited on 2023-03-07 02:58:04 by stolfi

PROG = test_frgb_path

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-tests 

all: check

check:  clean run-tests

run-tests: ${PROG} Makefile
	${PROG}
	./show_results.sh 

        
clean::
	rm -f out/*.txt out/*.ppm out/*.png
