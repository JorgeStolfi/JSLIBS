# Last edited on 2016-12-26 22:20:13 by stolfilocal

PROG := fn1_graphs
  
TEST_LIB := libfgraph.a
TEST_LIB_DIR := ../..

OTHER_I_FLAGS := \
  -IJSLIBS/libjs \
  -IJSLIBS/libps \
  -IJSLIBS/libflt \
  -IJSLIBS/libia \
  -IJSLIBS/libaa

OTHER_LIBS := \
  JSLIBS/libaa/libaa.a \
  JSLIBS/libia/libia.a \
  JSLIBS/libflt/libflt.a \
  JSLIBS/libps/libps.a \
  JSLIBS/libjs/libjs.a

all: run
  
include GENERIC-LIB-TEST.make  

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

PSVIEW = evince
  
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
