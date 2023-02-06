# Last edited on 2023-02-04 00:08:40 by stolfi

PROG := test_povray_camera
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS:= \
  libimg.a \
  libgeo.a \
  libjs.a

all: check

check:  build run_a

# OTHER_C_FLAGS:=-pg
# OTHER_LD_FLAGS:=-pg

include ${HOME}/programs/c/GENERIC-LIB-TEST.make

run_a: ${PROG} run-povray.sh
	./${PROG} 
	./run-povray.sh
	display out/crosses.pgm
