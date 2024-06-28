# Last edited on 2024-06-22 21:01:32 by stolfi

PROG := bbopt2

TEST_LIB := libbbopt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libaagraph.a \
  libeps.a \
  libaa.a \
  libia.a \
  libflt.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

FUNCTIONS := \
  f3_ia
  
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

OUTNAME := out/bb2_${FUNCTION}
  
EPSFILE := ${OUTNAME}.eps

PSVIEW := evince
  
single: ${EPSFILE} show-ps

clean-single: 
	/bin/rm -f ${OUTNAME}*.eps

show-ps: ${EPSFILE}
	${PSVIEW} ${EPSFILE}
  
# Generating the plots:

${EPSFILE}: ${PROG}
	${PROG} ${FUNCTION}
  
endif
# ${FUNCTION}
######################################################################
