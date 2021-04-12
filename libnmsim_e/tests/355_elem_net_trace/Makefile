# Last edited on 2020-12-17 09:31:27 by jstolfi

PROG := nmsim_test_355_elem_net_trace
 
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
	mkdir -p out
	nice ${PROG}

do-clean:
	rm -f out/sim_*.txt
