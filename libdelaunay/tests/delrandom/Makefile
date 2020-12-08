# Last edited on 2011-12-25 02:17:53 by stolfi

PROG := delrandom

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
# PSVIEW := gv
PSVIEW := okular

do-test: ${PSFILE} ${EPSFILE}

clean::
	rm -f ${PSFILE} ${EPSFILE}
  

${PSFILE}: ${PROG}
	${PROG} 80 uniform ps
	${PSVIEW} ${PSFILE}
  
${EPSFILE}: ${PROG}
	${PROG} 12 normal eps
	${PSVIEW} ${EPSFILE}

