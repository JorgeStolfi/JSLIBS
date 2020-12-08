# Last edited on 2020-12-05 07:52:07 by jstolfi

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
	nice ${PROG}
	for ff in out/sim_*_elem_n[0-9]*.txt ; do \
	  ./plot_potentials.sh $$ff ;\
	done

do-clean:
	rm -f out/sim_*.txt
