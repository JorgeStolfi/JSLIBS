# Last edited on 2023-02-18 09:50:55 by stolfi

PROG := fn1_zeros
 
TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libaagraph.a \
  libaafuncs.a \
  libaa.a \
  libia.a \
  libflt.a \
  libeps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: run run-single show-ps clean-single
.SUFFIXES: 

# all: clean run
all: run

FUNCTIONS := \
  f5

NOFUNCTIONS := \
  f1 f2 f3 f4 \
  gsin3

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

OUT_PREFIX := out/zf1_${FUNCTION}

MAIN_EPS = ${OUT_PREFIX}_aa_p.eps

run-single: ${MAIN_EPS} show-ps

clean-single: 
	/bin/rm -f ${PROG}.o ${PROG} ${OUT_PREFIX}_*.eps
        
PLOTSTEPS := 400  

${MAIN_EPS}: ${PROG}
	${PROG} ${FUNCTION} ${PLOTSTEPS}
  
# PSVIEW = ghostview
PSVIEW := evince

show-ps: ${MAIN_EPS}  
	for arith in ia id ar aa ; do \
	  ${PSVIEW} ${OUT_PREFIX}_$$arith_*.eps ; \     
	done
        
endif
# ${FUNCTION}
######################################################################
