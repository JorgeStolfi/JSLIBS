# Last edited on 2020-10-11 15:23:06 by jstolfi

PROG := bbopt1

TEST_LIB := libbbopt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libfgraph.a \
  libps.a \
  libaa.a \
  libia.a \
  libflt.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

FUNCTIONS := \
  f1_ia \
  f2_ia
  
NONFUNCTIONS :=

.PHONY:: do-test

all: check

check:  clean do-test

do-test: ${PROG}
	@for fn in ${FUNCTIONS} ; do \
          make FUNCTION=$$fn single ; \
        done

clean::
	@for fn in ${FUNCTIONS} ; do \
          make FUNCTION=$$fn clean-single ; \
        done
  
######################################################################
# For recursive "make" of single function:
# Caller must define ${FUNCTION}

FUNCTION := FUNCTION.IS.UNDEFINED
ifneq "/${FUNCTION}" "/FUNCTION.IS.UNDEFINED"

OUTNAME := bb1-out-${FUNCTION}
  
PSFILE := ${OUTNAME}-doc.ps
EPSFILE := ${OUTNAME}-000001.eps

PSFILES := \
  ${PSFILE} \
  ${EPSFILE}

MAINPS := ${EPSFILE}
  
PSVIEW := okular
  
single: ${MAINPS} show-ps

clean-single: 
	/bin/rm -f ${PSFILES}

show-ps: ${MAINPS}
	${PSVIEW} ${PSFILE}
  
# Generating the plots:

${MAINPS}: ${PROG}
	${PROG} ${FUNCTION}
  
endif
# ${FUNCTION}
######################################################################
