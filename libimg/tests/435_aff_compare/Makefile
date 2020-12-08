# Last edited on 2020-11-06 15:22:47 by jstolfi

PROG = test_aff_compare

# Test runs
# Format: {PREFIX}/{NS}/{IMG1}/{A1}/{R1}/{IMG2}/{A2}/{R2}/{DEFORM}
# Where {PREFIX} is the output file name prefix, {IMG1} and {IMG2}
# are image names without dir and ext, {A1} and {A2} each
# are six numbers separated by ":", {R1} and {R2} are extra mag factors 
# for {A1} and {A2}, and {DEFORM}
# is the variable deformation to apply: "xytrans", "rotmag", 
# "xystretch", or "xyshear".

NS := 25

TESTS := \
  rotoff_r15/1.00/scan-01/1:0:0:1:42.15:33.39/15/scan-03/+0.50000:+0.86603:-0.86603:+0.50000:76.85:31.61/15/rotmag 

MORE_TESTS := \
  ident_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-01/1:0:0:1:64.00:64.00/10/xytrans \
  ident_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-01/1:0:0:1:64.00:64.00/10/rotmag \
  ident_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-01/1:0:0:1:64.00:64.00/10/xystretch \
  ident_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-01/1:0:0:1:64.00:64.00/10/xyshear \
   \
  rescan_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-02/1:0:0:1:63.20:65.20/10/xytrans \
   \
  rotated_r10/1.00/scan-01/1:0:0:1:64.00:64.00/9.5/scan-03/+0.52992:+0.84805:-0.84805:+0.52992:62.35:64.76/11.0/rotmag \
  rotated_r05/1.00/scan-01/1:0:0:1:64.00:64.00/4.5/scan-03/+0.52992:+0.84805:-0.84805:+0.52992:62.35:64.76/5.3/rotmag \
  rotated_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-03/+0.50000:+0.86603:-0.86603:+0.50000:62.35:64.76/10/xytrans \
   \
  rotoff_r15/1.00/scan-01/1:0:0:1:42.15:33.39/15/scan-03/+0.50000:+0.86603:-0.86603:+0.50000:76.85:31.61/14.5/rotmag \
  rotoff_r15/1.00/scan-01/1:0:0:1:42.15:33.39/15/scan-03/+0.50000:+0.86603:-0.86603:+0.50000:76.85:31.61/14.5/xytrans \
   \
  magnified_r5/1.00/scan-01/1:0:0:1:64.00:64.00/5/scan-04/1:0:0:1:62.00:62.00/10/xytrans \
  magnified_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-04/1:0:0:1:62.00:62.00/20/xytrans \
   \
  reduced_r5/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-05/1:0:0:1:64.50:64.20/5/xytrans \
  reduced_r10/0.50/scan-01/1:0:0:1:64.00:64.00/20/scan-05/1:0:0:1:64.50:64.20/10/xytrans \
   \
  persp_r5/2.00/scan-01/1:0:0:1:64.00:64.00/5/scan-06/1:0:0:1:65.60:63.20/5/xytrans \
   \
  persp_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-06/+1.1338:-0.4113:-0.0161:+0.7384:65.17:64.68/10.10/rotmag \
  persp_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-06/+1.1338:-0.4113:-0.0161:+0.7384:65.60:64.26/10.10/xystretch \
  persp_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-06/+1.1338:-0.4113:-0.0161:+0.7384:65.60:64.26/10.10/xyshear \
  persp_r10/1.00/scan-01/1:0:0:1:64.00:64.00/10/scan-06/+1.1338:-0.4113:-0.0161:+0.7384:65.17:64.68/10.10/xytrans \

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
UVMAX :=  ${word 2, ${TEST_ARGS}}

IMG1 :=   ${word 3, ${TEST_ARGS}}
A1 :=     ${subst :, ,${word 4, ${TEST_ARGS}}}
R1 :=     ${word 5, ${TEST_ARGS}}

IMG2 :=   ${word 6, ${TEST_ARGS}}
A2 :=     ${subst :, ,${word 7, ${TEST_ARGS}}}
R2 :=     ${word 8, ${TEST_ARGS}}

DEFORM := ${word 9, ${TEST_ARGS}}
               
single-test: ${PROG}
	${PROG} ${PREFIX} ${DEFORM} ${NS} ${UVMAX} \
          ${IMG1} ${A1} ${R1} \
          ${IMG2} ${A2} ${R2}
	plot-mismatch.sh ${PREFIX} ${DEFORM}

endif
# End ${TEST} section
########################################################################
