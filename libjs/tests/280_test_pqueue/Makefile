# Last edited on 2023-02-04 07:06:00 by stolfi

PROG := test_pqueue

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libeps.a \
  libgeo.a \
  libjs.a
  
TESTS := \
  000020 \
  000200 \
  002000
 
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
NITEMS :=  ${word 1, ${TEST_FIELDS}}
  
PREFIX := out/${subst :,-,${TEST}}
PSVIEW := evince

single:${PROG}
	${PROG} ${NITEMS}

# End ${TEST} section
########################################################################
endif

clean::
	rm -f out/*.ps out/*.eps 
