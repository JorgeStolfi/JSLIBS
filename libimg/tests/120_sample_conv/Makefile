# Last edited on 2010-08-14 23:11:15 by stolfi

PROG = test_sample_conv

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminu.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-tests

all: check

check:  clean do-tests

do-tests: 
	mkdir -p out
	${MAKE} MAXVAL=255 ISMASK=0 do-test
	${MAKE} MAXVAL=256 ISMASK=0 do-test
	${MAKE} MAXVAL=255 ISMASK=1 do-test
	${MAKE} MAXVAL=256 ISMASK=1 do-test
        
ifneq "/${MAXVAL}" "/"   

OUTFILE := out/test-m${MAXVAL}-i${ISMASK}.txt

do-test: ${PROG}
	- ${PROG} -maxval ${MAXVAL} -isMask ${ISMASK} -range 2.0 3.0 > ${OUTFILE} 2>&1
	cat ${OUTFILE}

clean::
	/bin/rm -fv out/test-*.ppm

endif
