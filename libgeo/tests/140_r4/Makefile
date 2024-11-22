# Last edited on 2024-11-11 11:24:54 by stolfi

PROG = test_r4
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test check-cover

all: check check-cover

check:  do-test

do-test: ${PROG}
	${PROG}

CHECK_COVER_SCRIPT := ../../../compare_h_with_test_c.sh
check-cover: ${PROG}.c ${TEST_LIB_DIR}/${TEST_LIB}
	${CHECK_COVER_SCRIPT} ${PROG} r4 r4x4 
