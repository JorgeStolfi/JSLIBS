# Last edited on 2020-12-13 22:09:31 by jstolfi

PROG := nmsim_test_360_elem_net_group_stats
 
TEST_LIB := libnmsim.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	mkdir -p out
	${PROG}
