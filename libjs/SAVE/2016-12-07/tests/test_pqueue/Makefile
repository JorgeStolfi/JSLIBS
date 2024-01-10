# Last edited on 2011-05-29 10:48:03 by stolfi

PROG := test_pqueue

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libps.a \
  libgeo.a \
  libjs.a
  
TESTS := \
  000020 \
  000200 \
  002000
 
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
NITEMS :=  ${word 1, ${TEST_FIELDS}}
  
PREFIX := out/${subst :,-,${TEST}}
PSVIEW := okular

single:${PROG}
	${PROG} ${NITEMS}

# End ${TEST} section
########################################################################
endif

clean::
	rm -f out/*.ps out/*.eps 
