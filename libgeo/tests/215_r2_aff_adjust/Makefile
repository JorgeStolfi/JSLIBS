# Last edited on 2020-10-16 03:46:18 by jstolfi

PROG = r2_aff_adjust_test

# Test runs
# Format: ${TAG}/${AM00}:${AM01}:${AM10}:${AM11}:${AD0}:${AD1}/${RM00}:${RM01}:${RM10}:${RM11}:${RD0}:${RD1}
# where ${AM00}...${AD1} is the target optimum, ${RM00}..${RD1} is the max deviation in 
# each element.

TESTS := \
  scaleonly/10:2:-2:10:20:30/5:0:0:5:0:0 \
  shearonly/10:2:-2:10:20:30/0:5:5:0:0:0 \
  matonly/10:2:-2:10:20:30/5:5:5:5:0:0 \
  everything/10:2:-2:10:20:30/5:5:5:5:5:5

ALL_TESTS := \
  disponly/10:2:-3:12:20:30/0:0:0:0:5:5 \
  scaleonly/10:2:-2:10:20:30/5:0:0:5:0:0 \
  shearonly/10:2:-2:10:20:30/0:5:5:0:0:0 \
  matonly/10:2:-2:10:20:30/5:5:5:5:0:0 \
  everything/10:2:-2:10:20:30/5:5:5:5:5:5

TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminn.a \
  libminu.a \
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
        
########################################################################
# Section for recursive "make".
# Client must define ${TEST}

ifneq "/${TEST}" "/"

TEST_ARGS := ${subst /, ,${TEST}}
TAG  :=  ${word 1, ${TEST_ARGS}}
AOPT :=  ${word 2, ${TEST_ARGS}}
RAD  :=  ${word 3, ${TEST_ARGS}}

NPLOT := 20

OUT_PREFIX := out/${TAG}
               
single-test: ${PROG}
	${PROG} \
          -outPrefix ${OUT_PREFIX} \
          -optimum ${subst :, ,${AOPT}} \
          -deviation ${subst :, ,${RAD}} \
          -nPlot ${NPLOT}
	for mm in quad ; do \
          plot-mismatch.sh ${OUT_PREFIX} $$mm ; \
        done
	

endif
# End ${TEST} section
########################################################################
