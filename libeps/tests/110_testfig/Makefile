

PROG := testfig

TEST_LIB := libeps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test show-epsf

all: check

check:  do-test

OUTDIR := out

EPSFILES := \
  ${OUTDIR}/fig_A.eps \
  ${OUTDIR}/fig_B.eps

# PSVIEW := ghostview
# PSVIEW := gv
PSVIEW := evince
        
do-test: ${EPSFILES} show-epsf
  
${EPSFILES} : ${PROG} ${OUTDIR}
	${PROG}

show-epsf: ${EPSFILES}
	for f in ${EPSFILES} ; do \
	  ${PSVIEW} $$f ; \
        done
        
${OUTDIR}:
	mkdir -p ${OUTDIR}

clean::
	-/bin/rm -f ${EPSFILES}
  
