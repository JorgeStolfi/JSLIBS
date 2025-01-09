# Last edited on 2025-01-08 01:23:27 by stolfi

PROG = test_integrate_recursive

TEST_LIB := libpstzin0.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a \
  libjsextra.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-tests

all: check

check:  clean do-tests

ZFUNCS := 01 02 03 04 05 11 

do-tests: ${PROG}
	mkdir -p out
	for nf in ${ZFUNCS} ; do \
	  ./run_test.sh $$nf ; \
        done
 
clean::
	/bin/rm -fv out/*-*.fni
