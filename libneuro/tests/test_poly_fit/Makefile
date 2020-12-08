# Last edited on 2020-10-11 15:32:43 by jstolfi

PROG = test_poly_fit
 
TEST_LIB := libneuro.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-single-test do-all-tests

all: check

check:  do-all-tests

# Trial indices:
BEZPOLYS :=  000 001 002  003 004 005  006 007 008 009   010 011 012 013 014
RANPOLYS :=  015 016 017 018
# TRIALS := ${RANPOLYS} ${BEZPOLYS}
TRIALS := ${RANPOLYS}

do-all-tests: 
	for tt in ${TRIALS}  ; do \
          ${MAKE} TRIAL="$$tt" do-single-test ; \
        done

######################################################################
# Recursive make - client must define ${TRIAL}
ifneq "/${TRIAL}" "/"

do-single-test: ${PROG} out
	rm -f out/*_t${TRIAL}*.txt
	${PROG} ${TRIAL}
	show_results.sh ${TRIAL}
endif
# end ${TRIAL} section
######################################################################

out: 
	mkdir -p out
