# Last edited on 2020-12-17 10:58:56 by jstolfi

PROG := nmsim_test_360_elem_net_sim
 
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
	rm -fv out/*.{png,txt}
	nice ${PROG}
	for ff in out/sim_*_elem_ne[0-9]*_trace.txt ; do \
	  nmsim_plot_potentials.sh $$ff `nmsim_extract_neuron_indices.sh $$ff` 495 605 -90 10;\
	done

do-clean:
	rm -f out/sim_*.txt
