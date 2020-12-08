# Last edited on 2019-05-28 15:14:02 by jstolfi

PROG := nmsim_test_003_write_read
 
TEST_LIB := libnmsim.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  do-test

do-test: ${PROG}
	-rm -f out/file*.txt
	${PROG}
	for f in out/file0_*.txt ; do \
	  prdiff -Bb $${f} $${f/0_/1_} ; \
	done
