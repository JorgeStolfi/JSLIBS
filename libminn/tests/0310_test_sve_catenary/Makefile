# Last edited on 2024-01-11 07:12:57 by stolfi

PROG = test_sve_catenary
 
TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
TESTS := \
  008:04.00
 
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
NK := ${word 1,${TESTARGS}}
WD := ${word 2,${TESTARGS}}

OUT_PREFIX := out/${NK}

OUT_INI_DAT_FILE := ${OUT_PREFIX}-ini.dat
OUT_FIN_DAT_FILE := ${OUT_PREFIX}-fin.dat

OUT_INI_FNI_FILE := ${OUT_PREFIX}-ini-plt.fni
OUT_FIN_FNI_FILE := ${OUT_PREFIX}-fin-plt.fni

single-test: ${PROG} 
	time ${PROG} ${NK} ${WD}
	show-curve.sh ${OUT_INI_DAT_FILE} &
	show-plot.sh ${OUT_INI_FNI_FILE} &
	show-curve.sh ${OUT_FIN_DAT_FILE} &
	show-plot.sh ${OUT_FIN_FNI_FILE} &
	xmessage "OK"


endif
# end ${TEST} make
########################################################################

clean::
	-/bin/rm -f *-ot.dat
