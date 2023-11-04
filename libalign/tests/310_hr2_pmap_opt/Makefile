# Last edited on 2023-10-15 18:12:15 by stolfi

PROG = hr2_pmap_opt_test

TEST_LIB := libalign.a
TEST_LIB_DIR := ../..

JS_LIBS :=  \
  libminn.a \
  libminu.a \
  libgeo.a \
  libimg.a \
  libjs.a

# Test runs
# Format: ${TAG}/${A00}:${A01}:${A02}:${A10}:...:${A22}/${R00}:${R01}:${R02}:${R10}:...:${R22}
# where ${A00}...${A22} is the target optimum, ${R00}..${R22} is the max deviation in 
# each element.

TESTS := \
  scaleonly/1:10:20:0:10:2:0:-2:10/0:0:0:0:5:0:0:0:5 \
  shearonly/1:10:20:0:10:2:0:-2:10/0:0:0:0:0:5:0:5:0 \
  matonly/1:10:20:0:10:2:0:-2:10/0:0:0:0:5:5:0:5:5 \
  affine/1:10:20:0:10:2:0:-2:10/0:5:5:0:5:5:0:5:5 \
  persp/1:10:20:0:10:2:0:-2:10/0:0:0:5:0:0:5:0:0

ALL_TESTS := \
  disponly/1:12:20:0:10:2:0:-3:12/0:5:5:0:0:0:0:0:0 \
  scaleonly/1:10:20:0:10:2:0:-2:10/0:0:0:0:5:0:0:0:5 \
  shearonly/1:10:20:0:10:2:0:-2:10/0:0:0:0:0:5:0:5:0 \
  matonly/1:10:20:0:10:2:0:-2:10/0:0:0:0:5:5:0:5:5 \
  affine/1:10:20:0:10:2:0:-2:10/0:5:5:0:5:5:0:5:5 \
  persp/1:10:20:0:10:2:0:-2:10/0:0:0:5:0:0:5:0:0

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
