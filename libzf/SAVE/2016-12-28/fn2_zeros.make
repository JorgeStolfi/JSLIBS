# Tests of 2D zero finding (level curve, implicit curve) with AA and IA
# Last edited on 2016-12-28 23:23:21 by stolfilocal

PROG := fn2_zeros
 
TEST_LIB := libzf.a
TEST_LIB_DIR := ../..

OTHER_I_FLAGS := \
  -IJSLIBS/libjs \
  -IJSLIBS/libps \
  -IJSLIBS/libflt \
  -IJSLIBS/libia \
  -IJSLIBS/libaa \
  -IJSLIBS/libfgraph

OTHER_LIBS := \
  JSLIBS/libfgraph/libfgraph.a \
  JSLIBS/libaa/libaa.a \
  JSLIBS/libia/libia.a \
  JSLIBS/libflt/libflt.a \
  JSLIBS/libps/libps.a \
  JSLIBS/libjs/libjs.a

include GENERIC-LIB-TEST.make

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

OUTNAME := out/z2-${FUNCTION}-${METHOD}-${ORDER}
  
MAINPS := ${OUTNAME}-doc.ps

EPSFILES := \
  ${OUTNAME}-ia.eps \
  ${OUTNAME}-aa.eps

PSFILES := \
  ${MAINPS} \
  ${EPSFILES}
  
single: ${MAINPS} show-ps

clean-single: 
	/bin/rm -f ${PROG}.o ${PROG} ${PSFILES}
  
# PSVIEW := gv
PSVIEW := evince
  
${MAINPS}: ${PROG}
	${PROG} ${FUNCTION} ${METHOD} ${ORDER}

show-ps: ${MAINPS}
	${PSVIEW} ${MAINPS}

endif
endif
endif
# ${METHOD}, ${FUNCTION}, ${ORDER}
######################################################################
