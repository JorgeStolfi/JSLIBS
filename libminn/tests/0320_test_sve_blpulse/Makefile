# Last edited on 2023-02-27 11:02:13 by stolfi

PROG = test_sve_blpulse
 
TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libfftw3.so
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
TESTS := \
  063
 
.PHONY:: all-tests single test

all: check

check:  all-tests

all-tests: ${PROG}
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt single-test ; \
        done
	
########################################################################
# Recursive make -- caller must define ${TEST}

TEST = TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

TESTARGS := ${subst :, ,${TEST}}
NS := ${word 1,${TESTARGS}}

OUT_PREFIX := out/${NS}

OUT_INI_AMP_FILE := ${OUT_PREFIX}-ini-a.dat
OUT_INI_PWS_FILE := ${OUT_PREFIX}-ini-p.dat
OUT_INI_PLT_FILE := ${OUT_PREFIX}-ini-plt.fni

OUT_FIN_AMP_FILE := ${OUT_PREFIX}-fin-a.dat
OUT_FIN_PWS_FILE := ${OUT_PREFIX}-fin-p.dat
OUT_FIN_PLT_FILE := ${OUT_PREFIX}-fin-plt.fni

single-test: ${PROG} 
	time ${PROG} ${NS}
	for ver in ini fin ; do \
          ff=${OUT_PREFIX}-$$ver ; \
          ( show-pulse.sh $$ff-a.dat & ) ; \
          ( show-spectrum.sh $$ff-p.dat & ) ; \
          ( show-plot.sh $$ff-plt.fni ) ; \
  	  xmessage "OK" ; \
        done

endif
# end ${TEST} make
########################################################################

clean::
	-/bin/rm -f *-ot.dat
