# Last edited on 2019-01-11 00:12:50 by jstolfi

PROG := nmsim_test_320_elem_synapse
 
TEST_LIB := libnmsim.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}
