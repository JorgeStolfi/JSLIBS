# Last edited on 2022-02-28 11:52:52 by stolfi

PROG = r2_align_quadopt_test

TEST_LIB := libalign.a
TEST_LIB_DIR := ../..

JS_LIBS :=  \
  libminn.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: all-tests single-test

all: check

check:  clean all-tests single-test

all-tests: ${PROG}
	for tt in ${TESTS} ; do \
          ${MAKE} TEST="$$tt" single-test ; \
        done
        
DATAFILE := out/f.dat

single-test: ${PROG}
	-${PROG}
	@if [[ -r ${DATAFILE} ]]; then \
          for df in out/f*.dat ; do \
            plot_mismatch.sh $${df} ; \
          done ; \
        fi

clean:: 
	rm -f out/f*.dat out/f*.png
