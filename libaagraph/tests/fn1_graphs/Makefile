# Last edited on 2023-03-29 19:25:06 by stolfi

PROG := fn1_graphs
  
TEST_LIB := libaagraph.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libeps.a \
  libaafuncs.a \
  libaa.a \
  libia.a \
  libflt.a \
  libgeo.a \
  libjs.a
  
all: check

check:  run
  
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

OUTNAME := out/pf-${FUNCTION}

PLOTSTEPS := 512

MAINEPS := ${OUTNAME}-ia.eps

OTHEREPS :=  \
  ${filter-out ${MAINEPS},${wildcard ${OUTNAME}-*.eps}}

EPSFILES := \
  ${MAINEPS} \
  ${OTHEREPS}

EPSVIEW = evince
  
run-single: ${MAINEPS} show-eps

clean-single: 
	/bin/rm -f ${PROG}.o ${PROG} ${OUTNAME}-*.eps
  
${MAINEPS}: ${PROG}
	${PROG} ${FUNCTION} ${PLOTSTEPS}

show-eps: ${MAINEPS}        
	${EPSVIEW} ${OUTNAME}-*.eps
        
endif
# ${FUNCTION}
######################################################################
