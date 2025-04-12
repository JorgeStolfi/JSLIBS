# Last edited on 2025-03-19 12:51:42 by stolfi

PROG = test_r2_opt

# Test runs
# Format: ${FUNC}:${BIAS}:${MONO}:${QUAD}:${NI}:${RADX}:${RADY}:${STPX}:${STPY}:${SCLX}:${SCLY}

TESTS := \
  imgmis/0/1/1/02/02/02/00.25/00.25/0/0

PLOT_TEST := \
  imgmis/0/1/0/02/02/02/00.00/00.00/0/0

ALL_TESTS := \
  optdst/0/1/0/02/04/04/01.00/01.00/0/0 \
  optdst/0/1/1/02/04/04/01.00/01.00/0/0 \
  indiff/1/1/1/02/04/04/01.00/01.00/0/0 \
  imgmis/0/1/0/02/02/02/00.25/00.25/0/0 \
  imgmis/0/1/0/02/20/20/01.00/01.00/4/4 

TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminn.a \
  libminu.a \
  libgeo.a \
  libgausol.a \
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

TEST_ARGS := ${subst /, ,${TEST}}
FUNC :=   ${word  1, ${TEST_ARGS}}
BIAS :=   ${word  2, ${TEST_ARGS}}
MONO :=   ${word  3, ${TEST_ARGS}}
QUAD :=   ${word  4, ${TEST_ARGS}} 
NI :=     ${word  5, ${TEST_ARGS}}
RADX :=   ${word  6, ${TEST_ARGS}}
RADY :=   ${word  7, ${TEST_ARGS}}
STPX :=   ${word  8, ${TEST_ARGS}}
STPY :=   ${word  9, ${TEST_ARGS}}
SCLX :=   ${word 10, ${TEST_ARGS}}
SCLY :=   ${word 11, ${TEST_ARGS}}
              
single-test: ${PROG}
	-rm -f out/*.{png,pgm,dat}
	-${PROG} \
            ${FUNC}  ${BIAS} ${MONO} ${QUAD} ${NI}  \
            ${RADX} ${RADY}  ${STPX} ${STPY}  ${SCLX} ${SCLY}
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
