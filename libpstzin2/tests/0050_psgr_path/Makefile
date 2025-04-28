# Last edited on 2025-04-25 06:41:36 by stolfi

PROG = test_psgr_path

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
	cat out/path-n*.txt

CHECK_COVER_SCRIPT := ../../../compare_h_with_test_c.sh
COVER_HFILES := \
  psgr_path
  
check-cover: ${PROG}.c ${TEST_LIB_DIR}/${TEST_LIB}
	${CHECK_COVER_SCRIPT} ${PROG} ${COVER_HFILES}
        
clean::
	/bin/rm -fv out/*.txt
