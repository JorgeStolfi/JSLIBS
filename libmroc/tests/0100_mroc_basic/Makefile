# Last edited on 2021-07-08 02:26:23 by jstolfi

PROG = test_mroc_basic
 
TEST_LIB := libmroc.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

# Initial guess for number of triangles:
NF_GUESS := 10000

# Output file prefix:
OUT_PREFIX := out/basic

all: 
	mkdir -p out
	rm -f out/*.stl 
	${PROG} 
	stmesh_view -format ascii -eps 0.005 -nfGuess ${NF_GUESS} ${OUT_PREFIX}_tet.stl
	stmesh_view -format ascii -eps 0.005 -nfGuess ${NF_GUESS} ${OUT_PREFIX}_iso.stl

clean::
	rm -f out/*.stl 
