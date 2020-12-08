

PROG := teststream

TEST_LIB := libps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test show-ps

all: check

check:  do-test

OUTDIR := out

EPSFILES := \
  ${OUTDIR}/000001.eps \
  ${OUTDIR}/000002.eps \
  ${OUTDIR}/pgA.eps \
  ${OUTDIR}/pgB.eps
  
PSFILES := \
  ${OUTDIR}/doc-P.ps \
  ${OUTDIR}/doc-L.ps

# PSVIEW := ghostview
# PSVIEW := gv
PSVIEW := evince
        
do-test: ${PSFILES} ${EPSFILES} show-ps
  
${PSFILES} ${EPSFILES} : ${PROG} ${OUTDIR}
	rm -f ${OUTDIR}/*.{ps,eps}
	${PROG}

show-ps: ${PSFILES} ${EPSFILES}
	for f in ${PSFILES} ${EPSFILES} ; do \
	  ${PSVIEW} $$f ; \
        done
        
${OUTDIR}:
	mkdir -p ${OUTDIR}

clean::
	-/bin/rm -f ${PSFILE} ${EPSFILE}
  
