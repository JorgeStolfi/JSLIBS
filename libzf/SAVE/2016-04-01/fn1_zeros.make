# Last edited on 2011-05-29 10:52:21 by stolfi

PROG := fn1_zeros
 
TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libaa.a \
  libia.a \
  libflt.a \
  libfgraph.a \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: run run-single show-ps clean-single
.SUFFIXES: 

# all: clean run
all: run

FUNCTIONS := \
  gsin3

NOFUNCTIONS := \
  f1 f2 f3 f4 \
  f5

run: ${PROG}
	@for fn in ${FUNCTIONS} ; do \
          make FUNCTION=$$fn run-single ; \
        done

clean::
	@for fn in ${FUNCTIONS} ; do \
          make FUNCTION=$$fn clean-single ; \
        done
 
  
######################################################################
# For recursive "make" of single function:
# Caller must define ${FUNCTION}

FUNCTION = FUNCTION.IS.UNDEFINED
ifneq "/${FUNCTION}" "/FUNCTION.IS.UNDEFINED"

OUTNAME := z1-${FUNCTION}

ROOTSPS = ${OUTNAME}-s.ps
TRACEPS = ${OUTNAME}-p.ps

MAINPS = ${ROOTSPS}

PLOTSTEPS := 512
  
EPSFILES = \
  ${OUTNAME}-ia.eps \
  ${OUTNAME}-aa.eps

PSFILES = \
  ${MAINPS} \
  ${EPSFILES}
 
run-single: ${MAINPS} show-ps

clean-single: 
	/bin/rm -f ${PROG}.o ${PROG} ${PSFILES}
  
# PSVIEW = ghostview
PSVIEW := okular

${MAINPS}: ${PROG}
	${PROG} ${FUNCTION} ${PLOTSTEPS}

show-ps: ${MAINPS}        
	${PSVIEW} ${ROOTSPS}
	${PSVIEW} ${TRACEPS}
        
endif
# ${FUNCTION}
######################################################################
