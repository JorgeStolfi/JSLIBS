# Last edited on 2020-12-16 00:30:55 by jstolfi

PROG := nmsim_test_420_elem_neuron_trace_stats
 
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
