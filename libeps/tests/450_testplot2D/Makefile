# Last edited on 2021-04-13 23:23:22 by jstolfi

PROG := testplot2D
  
TEST_LIB := libpsextra.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libeps.a \
  libimg.a \
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
	${PROG}

show-eps: ${EPSFILES}        
	for f in  ${EPSFILES} ; do \
	  ${PSVIEW} $$f ; \
        done

clean::
	-/bin/rm -f ${EPSFILES}
