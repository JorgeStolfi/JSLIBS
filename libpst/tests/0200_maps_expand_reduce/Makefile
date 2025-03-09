# Last edited on 2025-02-27 23:03:21 by stolfi

PROG = test_maps_shrink_expand

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

# ZFUNCS := 01-ramp10 02-ramp01 04-parabo 05-spdom1 10-hcliff
ZFUNCS := 03-ramp11

do-tests: ${PROG}
	mkdir -p out
	for fp in ${ZFUNCS} ; do \
	  ./run_test.sh $$fp ; \
        done
 
clean::
	/bin/rm -fv out/*-*.{fni.png}
