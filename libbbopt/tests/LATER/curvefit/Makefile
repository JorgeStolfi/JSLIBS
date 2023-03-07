# Last edited on 2023-02-04 07:06:20 by stolfi

PROG := curvefit

TEST_LIB := libbbopt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libeps.a \
  libaa.a \
  libia.a \
  libjs.a \
  libflt.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test
  
OUTNAME := cft

PSFILE := ${OUTNAME}.ps
EPSFILE := ${OUTNAME}-000001.eps

PSFILES := \
  ${PSFILE} \
  ${EPSFILE}

MAINPS := ${EPSFILE}
  
PSVIEW := evince
  
do-test: ${MAINPS} show-ps

clean:: 
	/bin/rm -f ${PSFILES}

show-ps: ${MAINPS}
	${PSVIEW} ${PSFILE}
  
${MAINPS}: ${PROG}
	${PROG}
