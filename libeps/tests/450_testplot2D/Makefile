# Last edited on 2023-02-27 20:15:52 by stolfi

PROG := testplot2D
  
TEST_LIB := libeps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libeps.a \
  libimg.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run show-ps

all: check

check:  run

OUTDIR := out

EPSFILES := \
  ${OUTDIR}/plot2D_tri.eps \
  ${OUTDIR}/plot2D_qua.eps

PSVIEW = evince
        
run: ${PROG} show-eps
	-rm -rf ${EPSFILES}

${EPSFILES}: ${PROG} Makefile 
	${PROG}

show-eps: ${EPSFILES}        
	for f in  ${EPSFILES} ; do \
	  ${PSVIEW} $$f ; \
        done

clean::
	-/bin/rm -f ${EPSFILES}
