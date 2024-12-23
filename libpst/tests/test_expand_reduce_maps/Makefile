# Last edited on 2024-12-23 07:13:46 by stolfi

PROG = test_expand_reduce_maps

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
	  ./run-test.sh $$nf ; \
        done
 
clean::
	/bin/rm -fv out/*-*.fni
