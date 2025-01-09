# Last edited on 2025-01-05 10:53:23 by stolfi

PROG = test_pst_img_graph

TEST_LIB := libpstzin2.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
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
