# Last edited on 2023-02-20 06:23:37 by stolfi

PROG := delrandom

TEST_LIB := libdelaunay.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libquad.a \
  libeps.a \
  libgeo.a \
  libjs.a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  
	${MAKE} NSITES=60 DISTR=uniform do-test
	${MAKE} NSITES=60 DISTR=normal do-test
  
# PSVIEW := gv
PSVIEW := evince

# ######################################################################
# The ${NSITES} ${DISTR} section
ifneq "/${NSITES}" "/"
ifneq "/${DISTR}" "/"

EPSFILE = ${shell printf "out/delrandom_%05d_%s.eps" ${NSITES} "${DISTR}"}

do-test: ${EPSFILE}

clean::
	rm -f ${EPSFILE}

${EPSFILE}: ${PROG}
	${PROG} ${NSITES} ${DISTR}
	${PSVIEW} ${EPSFILE}

endif
endif
# End ${NSITES} ${DISTR} section
# ######################################################################
