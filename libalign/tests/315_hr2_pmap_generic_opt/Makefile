# Last edited on 2024-11-08 11:17:42 by stolfi

PROG = test_hr2_pmap_generic_opt

TEST_LIB := libalign.a
TEST_LIB_DIR := ../..

JS_LIBS :=  \
  libminn.a \
  libminu.a \
  libgeo.a \
  libimg.a \
  libjs.a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test
 
all: check

check:  do-test

OUT_PREFIX := out/test

do-test: ${PROG}
	${PROG}
	for datafile in ${OUT_PREFIX}-*-1D-plot.txt ; do \
          plot_opt_mismatch_1D.sh $$datafile ; \
        done
