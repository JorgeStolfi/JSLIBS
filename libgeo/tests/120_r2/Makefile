# Last edited on 2024-11-11 11:20:36 by stolfi

PROG = test_r2
 
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
	${CHECK_COVER_SCRIPT} ${PROG} r2 r2_extra r2_bezier r2x2
