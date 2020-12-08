# Last edited on 2020-11-06 13:04:43 by jstolfi

PROG = test_aff_extract

# Test runs
# Format: {PREFIX}/{IMG}/{A}/{AMAG}
# Where {PREFIX} is the output file name prefix, {IMG}
# is an image name without dir and ext, {A} is six numbers 
# separated by ":", and {AMAG} is an extra mag factor 
# for {A}.

TESTS := \
  rotoff_r10/scan-03/+0.52992:+0.84805:-0.84805:+0.52992:72.00:22.00/11.0 

MORE_TESTS := \
  ident_r10/scan-01/1:0:0:1:64.00:64.00/10 \
  ident_r20/scan-01/1:0:0:1:64.00:64.00/20 \
  ident_r05/scan-01/1:0:0:1:64.00:64.00/5 \
   \
  offset_r10/scan-01/1:0:0:1:47.00:43.00/10 \
  rotoff_r10/scan-03/+0.52992:+0.84805:-0.84805:+0.52992:72.00:22.00/11.0  \
   \
  rotated_r05/scan-03/+0.52992:+0.84805:-0.84805:+0.52992:62.35:64.76/5.3 \
  rotated_r10/scan-03/+0.52992:+0.84805:-0.84805:+0.52992:62.35:64.76/11.0 \
   \
  magnified_r05/scan-04/1:0:0:1:62.00:62.00/10 \
  magnified_r10/scan-04/1:0:0:1:62.00:62.00/20 \
   \
  reduced_r05/scan-05/1:0:0:1:64.50:64.20/5 \
  reduced_r10/scan-05/1:0:0:1:64.50:64.20/10 \
   \
  persp_r5/scan-06/1:0:0:1:65.60:63.20/5 \
   \
  persp_r10/scan-06/+1.1338:-0.4113:-0.0161:+0.7384:65.60:64.26/10.10

ALL_TESTS :=

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminn.a \
  libminu.a \
  libgeo.a \
  libjs.a 
  
OTHER_LIBS := \
  /usr/lib/libnetpbm.so \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng16.so \
  /usr/lib/x86_64-linux-gnu/libz.so

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
PREFIX := ${word 1, ${TEST_ARGS}}

IMG :=    ${word 2, ${TEST_ARGS}}
A :=      ${subst :, ,${word 3, ${TEST_ARGS}}}
AMAG :=   ${word 4, ${TEST_ARGS}}

single-test: ${PROG}
	${PROG} ${PREFIX} ${IMG} ${A} ${AMAG}
	display -filter Box -resize '400%' out/${PREFIX}_${IMG}.png

endif
# End ${TEST} section
########################################################################
