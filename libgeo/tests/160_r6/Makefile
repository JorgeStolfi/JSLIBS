# Last edited on 2024-11-11 11:31:38 by stolfi

PROG = test_r6
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test check-cover

all: check check-cover

check:  do-test

build-test: ${PROG}

do-test: ${PROG}
	${PROG}

CHECK_COVER_SCRIPT := ../../../compare_h_with_test_c.sh
check-cover: ${PROG}.c ${TEST_LIB_DIR}/${TEST_LIB}
	${CHECK_COVER_SCRIPT} ${PROG} r6 r6x6 
