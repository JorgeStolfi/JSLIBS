
PROG := testshade

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

EPSFILE24 := ${OUTDIR}/fig_002_004.eps

# PSVIEW := ghostview
# PSVIEW := gv
PSVIEW := evince
        
run: ${EPSFILE24} show-epsf
  
${EPSFILE24} : ${PROG}
	${PROG}

show-epsf: ${EPSFILE24}        
	${PSVIEW} ${OUTDIR}/*.eps

clean::
	-/bin/rm -f ${OUTDIR}/*.eps
