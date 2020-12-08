# Last edited on 2017-10-26 22:31:10 by stolfilocal

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_multifok_op

JS_LIBS := \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  ${PROG}
	mkdir -p out
	./${PROG} >& out/log.txt
