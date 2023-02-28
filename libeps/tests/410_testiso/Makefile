

PROG := testiso

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

EPSFILE := ${OUTDIR}/fig_003_003.eps

# PSVIEW := ghostview
# PSVIEW := gv
PSVIEW := evince
        
run: ${EPSFILE} show-epsf
  
${EPSFILE} : ${PROG}
	${PROG}

show-epsf: ${EPSFILE}        
	${PSVIEW} ${EPSFILE}

clean::
	-/bin/rm -f ${OUTDIR}/*.eps
