# Last edited on 2020-12-05 07:51:25 by jstolfi

PROG := nmsim_test_330_elem_neuron_trace
 
TEST_LIB := libnmsim_e.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libnmsim.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
  ${NN} ${NSTEPS}
all: check

check:  do-test

do-test: do-clean ${PROG}
	nice ${PROG}

do-clean:
	rm -f out/sim_*.txt
