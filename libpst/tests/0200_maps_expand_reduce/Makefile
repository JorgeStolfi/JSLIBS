# Last edited on 2024-12-23 09:16:29 by stolfi

PROG = test_maps_expand_reduce

TEST_LIB := libpst.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a \
  libjsextra.a

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
