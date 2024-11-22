# Last edited on 2023-02-18 11:17:09 by stolfi

PROG := test_mst

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libeps.a \
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
OUT_PREFIX := ${OUTDIR}/mst_${subst :,_,${TEST}}

# PSVIEW := gv
PSVIEW := atril

MAIN_EPS = ${OUT_PREFIX}.eps

single: ${MAIN_EPS}

${MAIN_EPS}: ${PROG} ${OUTDIR}
	${PROG} ${NSITES} ${DISTR} ${ACOST}
	${PSVIEW} ${MAIN_EPS}
        
${OUTDIR}:
	mkdir -p ${OUTDIR}

# End ${TEST} section
########################################################################
endif

clean::
	rm -f ${OUTDIR}/*.ps ${OUTDIR}/*.eps 
