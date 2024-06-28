# Tests of 2D zero finding (level curve, implicit curve) with AA and IA
# Last edited on 2024-06-22 20:55:19 by stolfi

PROG := fn2_zeros
 
TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libaagraph.a \
  libaafuncs.a \
  libeps.a \
  libaa.a \
  libia.a \
  libflt.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

FUNCTIONS := \
  f0 f1 f2 f3 f4 \
  f5 f6 f7 f8 \
  fadd \
  fdiv \
  fiamany \
  fmax \
  fmax2 \
  fmul \
  fmul2

NOFUNCTIONS :=

.PHONY:: run single show-ps clean-single

all: run

run: ${PROG}
	@for fn in ${FUNCTIONS} ; do \
          make FUNCTION=$$fn METHOD=quad ORDER=128 single ; \
          make FUNCTION=$$fn METHOD=grid ORDER=064 single ; \
        done

clean::
	@for fn in ${FUNCTIONS} ; do \
          echo make FUNCTION=$$fn METHOD=quad ORDER=128 clean-single ; \
          make FUNCTION=$$fn METHOD=grid ORDER=064 clean-single ; \
        done
  
######################################################################
# For recursive "make" of single function.
# Caller must define ${METHOD}, ${FUNCTION}, ${ORDER}.

METHOD := METHOD.IS.UNDEFINED
FUNCTION := FUNCTION.IS.UNDEFINED
ORDER := ORDER.IS.UNDEFINED

ifneq "/${METHOD}" "/METHOD.IS.UNDEFINED"
ifneq "/${FUNCTION}" "/FUNCTION.IS.UNDEFINED"
ifneq "/${ORDER}" "/ORDER.IS.UNDEFINED"

OUT_PREFIX := out/zf2_${FUNCTION}_${METHOD}_${ORDER}
  
MAIN_EPS := ${OUT_PREFIX}_quad_ia.eps
  
single: ${MAIN_EPS} show-ps

clean-single: 
	/bin/rm -f ${PROG}.o ${PROG} ${OUT_PREFIX}_*.eps
  
# PSVIEW := gv
PSVIEW := evince
  
${MAIN_EPS}: ${PROG}
	${PROG} ${FUNCTION} ${METHOD} ${ORDER}

show-ps: ${MAIN_EPS}
	for meth in quad grid ; do \
	  ${PSVIEW} ${OUT_PREFIX}_$$meth_${ORDER}_*.eps ; \
        done

endif
endif
endif
# ${METHOD}, ${FUNCTION}, ${ORDER}
######################################################################
