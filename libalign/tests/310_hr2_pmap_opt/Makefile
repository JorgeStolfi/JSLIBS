# Last edited on 2024-09-16 12:49:56 by stolfi

PROG = hr2_pmap_opt_test

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

# TYPE := TRANSLATION
TYPE := CONGRUENCE

OUT_PREFIX := out/test

do-test: ${PROG}
	${PROG} ${TYPE}
	for datafile in ${OUT_PREFIX}-*-${TYPE}-*-1D-plot.txt ; do \
          plot_mismatch_1D.sh $$datafile ; \
        done
