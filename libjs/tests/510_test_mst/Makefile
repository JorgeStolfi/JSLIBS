# Last edited on 2021-07-21 12:56:28 by jstolfi

PROG := test_mst

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libps.a \
  libgeo.a \
  libjs.a
  
TESTS := \
  0020:normal:edist \
  0020:normal:qdist \
  0200:normal:qdist \
  
 
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
ACOST :=   ${word 3, ${TEST_FIELDS}}
  
OUTDIR := out
PREFIX := ${OUTDIR}/${subst :,-,${TEST}}
EPSFILE = ${PREFIX}-000001.eps

# PSVIEW := gv
PSVIEW := atril

single: ${EPSFILE}

${EPSFILE}: ${PROG} ${OUTDIR}
	${PROG} ${NSITES} ${DISTR} ${ACOST} ${PREFIX}
	${PSVIEW} ${EPSFILE}
        
${OUTDIR}:
	mkdir -p ${OUTDIR}

# End ${TEST} section
########################################################################
endif

clean::
	rm -f ${OUTDIR}/*.ps ${OUTDIR}/*.eps 
