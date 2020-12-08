# Last edited on 2011-05-29 10:49:55 by stolfi

PROG := curvefit

TEST_LIB := libbbopt.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libps.a \
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
  
PSVIEW := okular
  
do-test: ${MAINPS} show-ps

clean:: 
	/bin/rm -f ${PSFILES}

show-ps: ${MAINPS}
	${PSVIEW} ${PSFILE}
  
${MAINPS}: ${PROG}
	${PROG}
