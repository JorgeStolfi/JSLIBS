# Last edited on 2024-11-06 10:27:01 by stolfi

PROG = test_aff_extract

# Test runs
# Format: {PREFIX}/{IMG}/{CTR}/{A}/{AMAG}
# Where {PREFIX} is the output file name prefix, {IMG}
# is an image name without dir and ext, {A} is six numbers 
# separated by ":", and {AMAG} is an extra mag factor 
# for {A}.

ALL_TESTS := \
  ident_r100/scan-01/64.00:64.00/1:0:0:1/1.00 \
  ident_r050/scan-01/64.00:64.00/1:0:0:1/0.50 \
  ident_r025/scan-01/64.00:64.00/1:0:0:1/0.25 \
   \
  offst_r010/scan-01/47.00:43.00/1:0:0:1/0.10  \
   \
  rotat_r025/scan-03/62.35:64.76/+0.52992:+0.84805:-0.84805:+0.52992/0.25 \
  rotat_r010/scan-03/62.35:64.76/+0.52992:+0.84805:-0.84805:+0.52992/0.10 \
  rotof_r010/scan-03/72.00:22.00/+0.52992:+0.84805:-0.84805:+0.52992/0.10 \
   \
  magnf_r010/scan-04/62.00:62.00/1:0:0:1/0.10 \
  magnf_r005/scan-04/62.00:62.00/1:0:0:1/0.05 \
   \
  reduc_r100/scan-05/64.50:64.20/1:0:0:1/1.00 \
  reduc_r200/scan-05/64.50:64.20/1:0:0:1/2.00 \
   \
  persp_r020/scan-06/65.60:63.20/1:0:0:1/0.20 \
   \
  persp_r099/scan-06/65.60:64.26/+1.1338:-0.4113:-0.0161:+0.7384/0.99

TESTS := ${ALL_TESTS}

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminn.a \
  libminu.a \
  libgeo.a \
  libjs.a 
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libnetpbm.so \
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
CTR :=    ${subst :, ,${word 3, ${TEST_ARGS}}}
A :=      ${subst :, ,${word 4, ${TEST_ARGS}}}
AMAG :=   ${word 5, ${TEST_ARGS}}

single-test: ${PROG}
	${PROG} ${PREFIX} ${IMG} ${CTR} ${A} ${AMAG}
	display -filter Box -resize '<400x' out/${PREFIX}_${IMG}.png

endif
# End ${TEST} section
########################################################################
