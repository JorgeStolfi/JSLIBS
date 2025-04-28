# Last edited on 2025-04-26 16:19:46 by stolfi

PROG = test_psgr_basic

TEST_LIB := libpstzin2.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
  libeps.a \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a \
  libjsextra.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test check-cover

all: check-cover check

check:  clean do-test

do-test: ${PROG}
	mkdir -p out
	${PROG}
	evince out/plot-lab0.eps
	evince out/plot-lab1.eps

CHECK_COVER_SCRIPT := ../../../compare_h_with_test_c.sh
COVER_HFILES := \
  psgr_types \
  psgr \
  psgr_copy \
  psgr_curl_image \
  psgr_plot \
  psgr_test \
  psgr_read \
  psgr_write
  
check-cover: ${PROG}.c ${TEST_LIB_DIR}/${TEST_LIB}
	${CHECK_COVER_SCRIPT} ${PROG} ${COVER_HFILES}
        
clean::
	/bin/rm -fv out/*-*.fni
