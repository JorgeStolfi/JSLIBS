# Last edited on 2025-01-04 21:58:46 by stolfi

PROG = test_virtual_gauge

TEST_LIB := libpst.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-tests

all: check

check:  clean do-tests

TESTS := 01 02
# TESTS := 01 02 03 04 05 11 

do-tests: ${PROG}
	mkdir -p out
	for nf in ${TESTS} ; do \
	  ./run_test.sh $$nf ; \
        done
 
clean::
	/bin/rm -fv out/*-*.fni
