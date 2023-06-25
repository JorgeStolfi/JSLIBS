# Last edited on 2023-03-27 09:45:38 by stolfi

PROG = test_minn

# Test runs
# Format: ${FUNC}:${N}

TESTS := \
  optnoq:1

PLOT_TEST := \
  optnoq:2

ALL_TESTS := \
  optdst:3 \
  optdst:3 \
  indiff:3 

TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminn.a \
  libminu.a \
  libgeo.a \
  libjs.a

.PHONY:: all-tests single-test

all: check

check:  clean all-tests single-test

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all-tests: ${PROG}
	for tt in ${TESTS} ; do \
          ${MAKE} TEST="$$tt" single-test ; \
        done
        
########################################################################
# Section for recursive "make".
# Client must define ${TEST}

ifneq "/${TEST}" "/"

TEST_ARGS := ${subst :, ,${TEST}}
FUNC :=   ${word  1, ${TEST_ARGS}}
N :=      ${word  2, ${TEST_ARGS}}
              
single-test: ${PROG}
	-rm -f out/*.{png,pgm,dat}
	-${PROG} ${FUNC} ${N}
	if [[ -s out/f2-00-00-ini.dat ]]; then \
	  for dd in `cd out && ls f2-*.dat` ; do plot_goal.sh $$dd ; done ; \
        fi
	if [[ -s out/image-00-00-000-ini.pgm ]]; then \
          display -title '%f' -filter Box -resize '3200%' out/image-*.pgm ; \
        fi

clean::
	rm -f out/f2*.*
        
endif
# End ${TEST} section
########################################################################
