# Last edited on 2024-11-06 11:16:13 by stolfi

PROG = test_lsq_robust
 
TEST_LIB := liblsq.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-single-test do-all-tests

all: check

check:  do-all-tests

# Trial indices:
BEZPOLYS :=  00:000  01:000 01:001  02:000 02:001 02:002  03:000 03:001 03:002 03:003   04:000 04:001 04:002 04:003 04:004
RANPOLYS :=  01:002  02:003 02:004  03:004 03:005         04:005 04:006                 05:006 05:007
BUGPOLYS :=                                               04:005 04:006                 05:006 05:007

# TRIALS := ${RANPOLYS} ${BEZPOLYS}
# TRIALS := ${RANPOLYS}
TRIALS := ${BUGPOLYS}

do-all-tests: 
	for tt in ${TRIALS}  ; do \
          ${MAKE} TRIAL="$$tt" do-single-test ; \
        done

######################################################################
# Recursive make - client must define ${TRIAL}
ifneq "/${TRIAL}" "/"

ARGS := ${subst :, ,${TRIAL}}
DEGREE := ${word 1, ${ARGS}}
IXPOLY  := ${word 2, ${ARGS}}

TAG := g${DEGREE}_ix${IXPOLY}

do-single-test: ${PROG} out
	rm -f out/*_${TAG}*.txt
	${PROG} ${DEGREE} ${IXPOLY}
	show_results.sh ${TAG}
endif
# end ${TRIAL} section
######################################################################

out: 
	mkdir -p out
