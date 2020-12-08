# Last edited on 2011-12-29 20:37:03 by stolfilocal

PROG := delgrid

TEST_LIB := libdelaunay.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libquad.a \
  libps.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  do-test

PSFILE = out/${PROG}-doc.ps
EPSFILE = out/${PROG}-000001.eps
PSVIEW := okular

do-test: ${PSFILE} ${EPSFILE}

clean::
	rm -f ${PSFILE} ${EPSFILE}
  

${PSFILE}: ${PROG}
	${PROG} 5 10 polar ps
	${PSVIEW} ${PSFILE}
  
${EPSFILE}: ${PROG}
	${PROG} 5 10 rect eps
	${PSVIEW} ${EPSFILE}

