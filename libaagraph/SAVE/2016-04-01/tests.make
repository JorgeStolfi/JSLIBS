# Last edited on 2011-05-29 10:48:29 by stolfi

PROG := fn1_graphs
  
TEST_LIB := libfgraph.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libps.a \
  libaa.a \
  libia.a \
  libflt.a \
  libjs.a
  
all: run
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make  

FUNCTIONS := \
  f4 f1 f2 f3 \

NONFUNCTIONS := \
  gexp \
  gasqrt \
  gsin3 \
  glog4 \
  ginv \
  giaboom \
  g1 g2 g3 g4 g5 g6 g7 g8 \
  g10 g11 g12 \
  g14 \
  gbadia \
  gbadmul \
  gdiv \
  gsqr \
  gsqrt
 
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

OUTNAME := pf-${FUNCTION}

PLOTSTEPS := 512

MAINPS := ${OUTNAME}-doc.ps

EPSFILES :=  \
  ${OUTNAME}-ia.eps \
  ${OUTNAME}-aa.eps \
  ${OUTNAME}-ar.eps

PSFILES := \
  ${MAINPS} \
  ${EPSFILES}

PSVIEW = okular
  
run-single: ${MAINPS} show-ps

clean-single: 
	/bin/rm -f ${PROG}.o ${PROG} ${PSFILES}
  
${MAINPS}: ${PROG}
	${PROG} ${FUNCTION} ${PLOTSTEPS}

show-ps: ${MAINPS}        
	${PSVIEW} ${MAINPS}
        
endif
# ${FUNCTION}
######################################################################
