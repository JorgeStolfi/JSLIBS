# Last edited on 2020-12-17 11:10:06 by jstolfi

PROG := nmsim_test_420_elem_neuron_trace_stats
 
TEST_LIB := libnmsim_e.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libnmsim.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  ${PROG}
	mkdir -p out
	rm -f out/*.{txt,png}
	nice ${PROG}
	
