# Last edited on 2019-01-02 22:54:29 by jstolfi

PROG := nmsim_test_005_firing_func
 
TEST_LIB := libnmsim.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}
	for ff in out/Phi_*.txt ; do \
	  ./plot_firing_func.sh $${ff%.*} ;\
	done
