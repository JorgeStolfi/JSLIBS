# Last edited on 2024-11-08 11:18:45 by stolfi

PROG = test_r2_align

TEST_LIB := libalign.a
TEST_LIB_DIR := ../..

JS_LIBS :=  \
  libminn.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: all-tests single-test

all: check

check:  clean single-test

DATAFILE := out/f005.dat

single-test: ${PROG}
	-rm -f out/f*.dat out/f*.png
	-${PROG} 
	@if [[ -r ${DATAFILE} ]]; then \
          for df in out/f*.dat ; do \
            plot_rel_dist_sqr.sh $${df} ; \
          done ; \
        fi

clean:: 
	rm -fv out/f*.dat out/f*.png
