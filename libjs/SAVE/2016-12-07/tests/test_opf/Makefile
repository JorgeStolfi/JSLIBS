# Last edited on 2011-05-29 10:48:19 by stolfi

PROG := test_opf

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libps.a \
  libgeo.a \
  libjs.a
  
TESTS := \
  0020:normal:fmax \
  0020:normal:fsum \
  0200:normal:fmax \
  
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: single

all: 
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
  
PREFIX := out/${subst :,-,${TEST}}
EPSFILE = ${PREFIX}-000001.eps
PSVIEW := okular

single: ${EPSFILE}

${EPSFILE}: ${PROG}
	${PROG} ${NSITES} ${DISTR} ${PCOST} ${PREFIX}
	${PSVIEW} ${EPSFILE}

# End ${TEST} section
########################################################################
endif

clean::
	rm -f out/*.ps out/*.eps 
