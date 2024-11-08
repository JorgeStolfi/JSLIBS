# Last edited on 2024-11-07 18:57:31 by stolfi

PROG = test_sve_img_match
 
TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
# Test run arguments: {TESTNUM}:{SHRINK}:{CTRX}:{CTRY}:{ANGLE}:{SCALE}
# 
# where {TESTNUM} is test index and the other parameters are 
# initial guesses. 

TESTS := \
  01:0001:0193:0122:023:1.33
  
NON_TESTS := \
  01:0010:0193:0122:023:1.33 \
  01:0001:0193:0122:023:1.33 \
   \
  01:0005:0173:0102:015:2.00

.PHONY:: all-tests single-test show-images

all: check

check:  all-tests

all-tests: ${PROG}
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt single-test show-images; \
        done

clean::
	-/bin/rm -f out/test-*.fni out/test-*.png 
	
########################################################################
# Recursive make -- caller must define ${TEST}

TEST = TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

TESTARGS := ${subst :, ,${TEST}}
TESTNUM :=  ${word 1, ${TESTARGS}}
SHRINK :=   ${word 2, ${TESTARGS}}
CTRX :=     ${word 3, ${TESTARGS}}
CTRY :=     ${word 4, ${TESTARGS}}
ANGLE :=    ${word 5, ${TESTARGS}}
SCALE :=    ${word 6, ${TESTARGS}}

IPREFIX := in/test-${TESTNUM}
IIMG := ${IPREFIX}-img
IOBJ := ${IPREFIX}-obj
IMSK := ${IPREFIX}-msk

OPREFIX := out/test-${TESTNUM}
OIMG_INI := ${OPREFIX}-ini-img
OOBJ_INI := ${OPREFIX}-ini-obj
OMSK_INI := ${OPREFIX}-ini-msk

OIMG_OPT := ${OPREFIX}-opt-img
OOBJ_OPT := ${OPREFIX}-opt-obj
OMSK_OPT := ${OPREFIX}-opt-msk

OPLOT := ${OPREFIX}-f2-plot

single-test: ${PROG} ${IIMG}.fni ${IOBJ}.fni ${IMSK}.fni 
	${PROG} \
            -img ${IIMG}.fni \
            -obj ${IOBJ}.fni \
            -msk ${IMSK}.fni \
            -shrink ${SHRINK} \
            -center ${CTRX} ${CTRY} \
            -angle ${ANGLE} \
            -scale ${SCALE} \
            -outPrefix ${OPREFIX}

show-images:
	cat ${OIMG_INI}.fni | fni_to_pnm -min 0 -max 1 | convert PGM:- ${OIMG_INI}.png
	cat ${OOBJ_INI}.fni | fni_to_pnm -min 0 -max 1 | convert PGM:- ${OOBJ_INI}.png
	cat ${OMSK_INI}.fni | fni_to_pnm -min 0 -max 1 | convert PGM:- ${OMSK_INI}.png
	cat ${OIMG_OPT}.fni | fni_to_pnm -min 0 -max 1 | convert PGM:- ${OIMG_OPT}.png
	cat ${OOBJ_OPT}.fni | fni_to_pnm -min 0 -max 1 | convert PGM:- ${OOBJ_OPT}.png
	cat ${OMSK_OPT}.fni | fni_to_pnm -min 0 -max 1 | convert PGM:- ${OMSK_OPT}.png
	display -title '%f' -filter box -resize '400%' \
	    ${OIMG_INI}.png ${OOBJ_INI}.png ${OMSK_INI}.png \
	    ${OIMG_OPT}.png ${OOBJ_OPT}.png ${OMSK_OPT}.png \
          & show-goalf-plot.sh ${OPLOT}.txt \
	  & ( xmessage "OK" ; echo "OK" 1>&2 )

${IIMG}.fni: ${IIMG}.png 
	convert ${IIMG}.png PGM:- | pnm_to_fni -min 0 -max 1 > ${IIMG}.fni

${IOBJ}.fni: ${IOBJ}.png 
	convert ${IOBJ}.png PGM:- | pnm_to_fni -min 0 -max 1 > ${IOBJ}.fni

${IMSK}.fni: ${IMSK}.png 
	convert ${IMSK}.png PGM:- | pnm_to_fni -min 0 -max 1 > ${IMSK}.fni

endif
# end ${TEST} make
########################################################################
