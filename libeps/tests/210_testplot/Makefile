

PROG := testplot
  
TEST_LIB := libeps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run show-epsf

all: check

check:  run

OUTDIR := out

EPSFILE := ${OUTDIR}/fig.eps

# PSVIEW = ghostview
# PSVIEW = gv
PSVIEW = evince
        
run: ${EPSFILE} show-epsf
  
${EPSFILE} : ${PROG}
	rm -f ${OUTDIR}/*.eps
	${PROG}

show-epsf: ${EPSFILE}        
	${PSVIEW} ${EPSFILE}

clean::
	-/bin/rm -f ${EPSFILE}
