# Last edited on 2024-11-08 11:16:37 by stolfi

PROG = test_r2_align_multiscale

TEST_LIB := libalign.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminn.a \
  libgeo.a \
  libimg.a \
  libjs.a

# Test runs
# Format: ${FUNC}:${BIAS}:${MONO}:${QUAD}:${NI}:${WX}:${WY}:${AX}:${AY}

TESTS := \
  points:0:0:0:02:04:04:01:01

ALL_TESTS := \
  images:0:0:0:02:20:20:01:01 \
  indiff:0:0:0:02:04:04:01:01 

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: all-tests single-test

all: check

check:  clean all-tests single-test

all-tests: ${PROG}
	for tt in ${TESTS} ; do \
          ${MAKE} TEST="$$tt" single-test ; \
        done
        
########################################################################
# Section for recursive "make".
# Client must define ${TEST}

ifneq "/${TEST}" "/"

TEST_ARGS := ${subst :, ,${TEST}}
FUNC := ${word 1, ${TEST_ARGS}}
BIAS := ${word 2, ${TEST_ARGS}}
MONO := ${word 3, ${TEST_ARGS}}
QUAD := ${word 4, ${TEST_ARGS}} 
NI :=   ${word 5, ${TEST_ARGS}}
WX :=   ${word 6, ${TEST_ARGS}}
WY :=   ${word 7, ${TEST_ARGS}}
AX :=   ${word 8, ${TEST_ARGS}}
AY :=   ${word 9, ${TEST_ARGS}}
               
single-test: ${PROG}
	-${PROG} ${FUNC} ${NI}  ${BIAS} ${MONO} ${QUAD}  ${WX} ${WY}  ${AX} ${AY} 
	-plot-mismatch.sh
	@if [[ -r out/img-00.]]; then \
          display -title '%f' -filter Point -resize '400%' `ls out/img-*.pgm | sort` ; \
        fi

endif
# End ${TEST} section
########################################################################
