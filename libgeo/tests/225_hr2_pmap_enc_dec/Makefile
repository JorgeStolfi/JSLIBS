# Last edited on 2024-09-17 19:41:39 by stolfi

PROG = hr2_pmap_enc_dec_test

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
