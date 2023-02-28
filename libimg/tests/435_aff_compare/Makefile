# Last edited on 2023-02-27 16:01:08 by stolfi

PROG = test_aff_compare

# Test runs
# Format: {PREFIX}/{NS}/{IMG1}/{A1}/{R1}/{IMG2}/{A2}/{R2}/{DEFORM}
# Where {PREFIX} is the output file name prefix, {IMG1} and {IMG2}
# are image names without dir and ext, {A1} and {A2} each
# are six numbers separated by ":", {R1} and {R2} are extra mag factors 
# for {A1} and {A2}, and {DEFORM}
# is the variable deformation to apply: "xytrans", "rotmag", 
# "xystretch", or "xyshear".

NS := 11

SOME_TESTS := \
  ident_r010_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-01/64.00:64.00/1:0:0:1/0.10/xytrans \
  ident_r010_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-01/64.00:64.00/1:0:0:1/0.10/xyshear \
  rotof_r011_009/1.00/scan-01/64.00:64.00/1:0:0:1/0.11/scan-03/62.35:64.76/+0.52992:+0.84805:-0.84805:+0.52992/0.09/rotmag \
  magnf_r020_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.20/scan-04/62.00:62.00/1:0:0:1/0.10/xytrans \
  persp_r010_009/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-06/65.60:64.26/+1.1338:-0.4113:-0.0161:+0.7384/0.09/xystretch

MORE_TESTS := \
  ident_r010_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-01/64.00:64.00/1:0:0:1/0.10/xytrans \
  ident_r010_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-01/64.00:64.00/1:0:0:1/0.10/rotmag \
  ident_r010_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-01/64.00:64.00/1:0:0:1/0.10/xystretch \
  ident_r010_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-01/64.00:64.00/1:0:0:1/0.10/xyshear \
   \
  rotof_r010_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-03/62.35:64.76/+0.50000:+0.86603:-0.86603:+0.50000/0.10/xytrans \
  rotof_r007_006/1.00/scan-01/42.15:33.39/1:0:0:1/0.07/scan-03/76.85:31.61/+0.50000:+0.86603:-0.86603:+0.50000/0.06/xytrans \
   \
  rotof_r011_009/1.00/scan-01/64.00:64.00/1:0:0:1/0.11/scan-03/62.35:64.76/+0.52992:+0.84805:-0.84805:+0.52992/0.09/rotmag \
  rotof_r022_018/1.00/scan-01/64.00:64.00/1:0:0:1/0.22/scan-03/62.35:64.76/+0.52992:+0.84805:-0.84805:+0.52992/0.18/rotmag \
  rotof_r007_008/1.00/scan-01/42.15:33.39/1:0:0:1/0.07/scan-03/76.85:31.61/+0.50000:+0.86603:-0.86603:+0.50000/0.08/rotmag \
  rotof_r007_007/1.00/scan-01/42.15:33.39/1:0:0:1/0.07/scan-03/76.85:31.61/+0.50000:+0.86603:-0.86603:+0.50000/0.07/rotmag \
   \
  magnf_r020_010/1.00/scan-01/64.00:64.00/1:0:0:1/0.20/scan-04/62.00:62.00/1:0:0:1/0.10/xytrans \
  magnf_r010_005/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-04/62.00:62.00/1:0:0:1/0.05/xytrans \
  magnf_r010_005/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-04/62.00:62.00/1:0:0:1/0.05/rotmag \
   \
  reduc_r010_020/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-05/64.50:64.20/1:0:0:1/0.20/xytrans \
  reduc_r005_010/0.50/scan-01/64.00:64.00/1:0:0:1/0.05/scan-05/64.50:64.20/1:0:0:1/0.10/xytrans \
  reduc_r005_010/0.50/scan-01/64.00:64.00/1:0:0:1/0.05/scan-05/64.50:64.20/1:0:0:1/0.10/rotmag \
   \
  persp_r010_009/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-06/65.17:64.68/+1.1338:-0.4113:-0.0161:+0.7384/0.09/rotmag \
  persp_r010_009/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-06/65.60:64.26/+1.1338:-0.4113:-0.0161:+0.7384/0.09/xystretch \
  persp_r010_009/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-06/65.60:64.26/+1.1338:-0.4113:-0.0161:+0.7384/0.09/xyshear \
  persp_r010_009/1.00/scan-01/64.00:64.00/1:0:0:1/0.10/scan-06/65.17:64.68/+1.1338:-0.4113:-0.0161:+0.7384/0.09/xytrans \

ALL_TESTS := ${SOME_TESTS} ${MORE_TESTS}

TESTS := ${SOME_TESTS}
# TESTS := ${ALL_TESTS}

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminn.a \
  libminu.a \
  libgeo.a \
  libjs.a 
  
OTHER_LIBS := \
  /usr/lib/libnetpbm.so.10 \
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
CTR1 :=   ${subst :, ,${word 4, ${TEST_ARGS}}}
A1 :=     ${subst :, ,${word 5, ${TEST_ARGS}}}
R1 :=     ${word 6, ${TEST_ARGS}}

IMG2 :=   ${word 7, ${TEST_ARGS}}
CTR2 :=   ${subst :, ,${word 8, ${TEST_ARGS}}}
A2 :=     ${subst :, ,${word 9, ${TEST_ARGS}}}
R2 :=     ${word 10, ${TEST_ARGS}}

DEFORM := ${word 11, ${TEST_ARGS}}
               
single-test: ${PROG}
	${PROG} ${PREFIX} ${DEFORM} ${NS} ${UVMAX} \
          ${IMG1} ${CTR1} ${A1} ${R1} \
          ${IMG2} ${CTR2} ${A2} ${R2}
	plot-mismatch.sh ${PREFIX} ${DEFORM}

endif
# End ${TEST} section
########################################################################
