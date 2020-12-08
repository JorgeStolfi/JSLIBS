# Last edited on 2019-04-09 19:45:01 by jstolfi

PROG := testshade

TEST_LIB := libps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run show-ps

all: check

check:  run

OUTDIR := out

EPSFILE := ${OUTDIR}/one.eps
PSFILE := ${OUTDIR}/doc.ps

# PSVIEW := ghostview
# PSVIEW := gv
PSVIEW := atril
        
run: ${PSFILE} ${EPSFILE} show-ps
  
${PSFILE} ${EPSFILE} : ${PROG}
	${PROG}

show-ps: ${PSFILE} ${EPSFILE}        
	${PSVIEW} ${PSFILE}
	${PSVIEW} ${EPSFILE}

clean::
	-/bin/rm -f ${PSFILE} ${EPSFILE}
