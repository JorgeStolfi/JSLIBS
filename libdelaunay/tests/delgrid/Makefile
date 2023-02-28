# Last edited on 2023-02-20 05:57:55 by stolfi

PROG := delgrid

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
	${MAKE} NI=5 NJ=10 GTYPE=polar do-test
	${MAKE} NI=7 NJ=7 GTYPE=rect do-test

PSVIEW := evince

# ######################################################################
# The ${NI} ${NJ} ${GTYPE} section
ifneq "/${NI}" "/"
ifneq "/${NJ}" "/"
ifneq "/${GTYPE}" "/"

NSITES = ${shell echo "${NI}*${NJ}" | bc -lq }
EPSFILE = ${shell printf "out/delgrid_%05d_%s.eps" ${NSITES} "${GTYPE}"}

do-test: ${EPSFILE}

clean::
	rm -f ${EPSFILE}

${EPSFILE}: ${PROG}
	${PROG} ${NI} ${NJ} ${GTYPE}
	${PSVIEW} ${EPSFILE}

endif
endif
endif
# End ${NI} ${NJ} ${GTYPE} section
# ######################################################################
