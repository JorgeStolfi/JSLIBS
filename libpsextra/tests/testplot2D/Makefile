# Last edited on 2011-06-06 17:51:01 by stolfi

PROG := testplot2D
  
TEST_LIB := libpsextra.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libps.a \
  libimg.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run show-ps

all: check

check:  run

OUTDIR := out

EPSFILES := \
  ${OUTDIR}/tri.eps \
  ${OUTDIR}/qua.eps

PSFILE := \
  ${OUTDIR}/doc.ps

# PSVIEW = ghostview
# PSVIEW = gv
PSVIEW = okular
        
run: ${PSFILE} ${EPSFILE} show-ps
  
${PSFILE}: ${PROG}
	-rm -rf ${PSFILE} ${EPSFILES}
	${PROG}

${EPSFILES}: ${PSFILE}

show-ps: ${PSFILE} ${EPSFILES}        
	${PSVIEW} ${PSFILE}
	for f in  ${EPSFILES} ; do \
	  ${PSVIEW} $$f ; \
        done

clean::
	-/bin/rm -f ${PSFILE} ${EPSFILES}
