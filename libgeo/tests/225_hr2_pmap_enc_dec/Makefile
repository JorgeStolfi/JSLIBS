# Last edited on 2024-11-08 11:27:07 by stolfi

PROG = test_hr2_pmap_enc_dec

TEST_LIB := libgeo.a
TEST_LIB_DIR := ../..

JS_LIBS :=  \
  libimg.a \
  libjs.a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test
 
all: check

check:  do-test

do-test: ${PROG}
	${PROG}
