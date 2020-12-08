# Last edited on 2018-06-30 04:40:35 by stolfilocal

PROG = test_expand_reduce_maps

TEST_LIB := libpst.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
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
	  ./run-test.sh $$nf ; \
        done
 
clean::
	/bin/rm -fv out/*-*.fni
