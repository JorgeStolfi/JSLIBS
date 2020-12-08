# Last edited on 2019-01-06 02:14:23 by jstolfi

PROG := nmsim_test_210_group_neuron
 
TEST_LIB := libnmsim.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}
