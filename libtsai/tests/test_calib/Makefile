# Last edited on 2023-02-03 23:21:11 by stolfi

PROG := test_calib
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS := \
  libtsai.a \
  libeps.a \
  libgeo.a \
  libjs.a

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg
 
include ../../../../GENERIC-LIB-TEST.make

all: check

check:  build run-tests

run-tests: t01

# Tests ----------------------------------------------------------------------

t01: ${PROG} in/t01
	-rm -rf out/t01/*
	./${PROG} \
	    optura \
	    in/t01 \
            out/t01 \
            2 1 0
	${PSVIEW} out/t01/*.eps

t02: ${PROG} in/t02
	-rm -rf out/t02/*
	./${PROG} \
	    svga \
	    in/t02 \
            out/t02 \
            2 1 1
	${PSVIEW} out/t02/*.eps
