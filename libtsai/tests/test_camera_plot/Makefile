# Last edited on 2023-02-04 07:05:34 by stolfi

PROG := test_camera_plot
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS:= \
  libgeo.a \
  libeps.a \
  libjs.a
  
all: check

check:  build run-tests

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

PSVIEW := evince

include ../../../../GENERIC-LIB-TEST.make

run-tests: ${PROG}
	-rm -rf out/*
	./${PROG}
	${PSVIEW} out/*.eps
