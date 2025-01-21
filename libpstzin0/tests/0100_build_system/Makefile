# Last edited on 2025-01-14 15:45:34 by stolfi

PROG = test_build_system

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

TESTS := 00-zeflat-N:T 01-ramp10-N:T 04-parabo-N:T 20-cpiece-N:T

do-tests: ${PROG}
	mkdir -p out
	for test in ${TESTS} ; do \
          mapName="$${test/:*/}"; full="$${test/*:/}" ; \
	  ./run_test.sh $${mapName} $${full} ; \
        done
 
clean::
	/bin/rm -fv out/*-*.fni
