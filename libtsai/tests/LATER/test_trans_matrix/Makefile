# Last edited on 2020-10-11 03:45:46 by jstolfi

PROG := test_trans_matrix
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS:= \
  liblmdif.a \
  libimg.a \
  libgeo.a \
  libjs.a

all: check

check:  build run-a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

run-%: out/test-%.pgm

out/test-%.pgm: ${PROG} ${wildcard data/test-%/*}
	rm -rf out/test-$*/
	./${PROG} \
	    data/test-$* \
            out/test-$*

