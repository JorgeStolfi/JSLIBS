# Last edited on 2023-02-18 13:17:07 by stolfi

PROG := test_opf

TEST_LIB := libgausol.a \
  libjs.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libeps.a \
  libgeo.a \
  libgausol.a \
  libjs.a
  
TESTS := \
  0020:normal:fmax \
  0020:normal:fsum \
  0200:normal:fmax \
  
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: single

all: check

check:  
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t single ; \
        done
        
ifneq "/${TEST}" "/"
########################################################################
# Section for single test, requires ${TEST}

TEST_FIELDS := ${subst :, ,${TEST}}
NSITES :=  ${word 1, ${TEST_FIELDS}}
DISTR :=   ${word 2, ${TEST_FIELDS}}
PCOST :=   ${word 3, ${TEST_FIELDS}}
  
OUT_PREFIX := out/opf_${subst :,_,${TEST}}
PSVIEW := evince

EPS_FILE := ${OUT_PREFIX}.eps

single: ${EPS_FILE}

${EPS_FILE}: ${PROG}
	${PROG} ${NSITES} ${DISTR} ${PCOST}
	${PSVIEW} ${EPS_FILE}

# End ${TEST} section
########################################################################
endif

clean::
	rm -f out/*.ps out/*.eps 
