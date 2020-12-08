# Last edited on 2011-05-15 20:22:46 by stolfi

PROG := test_camera_plot
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS:= \
  liblmdif.a \
  libgeo.a \
  libps.a \
  libjs.a
  
all: check

check:  build run-tests

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

PSVIEW := okular

include ../../../GENERIC-LIB-TEST.make

run-tests: ${PROG}
	-rm -rf out/*
	./${PROG}
	${PSVIEW} out/*.eps
