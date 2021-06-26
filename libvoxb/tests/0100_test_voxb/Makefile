# Last edited on 2021-06-11 16:22:53 by jstolfi

PROG = test_voxb
 
TEST_LIB := libvoxb.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libppv.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

TESTS := \
  0200:erolate

MORE_TESTS := \
  0200:objs \
  0200:tubes

all: check

check:  run-all-tests

run-all-tests:
	for test in ${TESTS} ; do \
          ${MAKE} TEST_ARGS="$${test}" run-single ; \
        done
        
######################################################################
# For recursive makes - caller must define ${TEST_ARGS}
ifneq "/${TEST_ARGS}" "/"

TEST_ARGS_SPLIT = ${subst :, ,${TEST_ARGS}}

SIZE :=   ${word 1, ${TEST_ARGS_SPLIT}}
OBJECT := ${word 2, ${TEST_ARGS_SPLIT}}

NF_GUESS := ${shell echo '80*${SIZE}*${SIZE}' | bc -lq}

TEST_NAME := tt_${SIZE}_${OBJECT}

OUT_PREFIX := out/${TEST_NAME}

TOM_FILE := ${OUT_PREFIX}.tom
STL_FILE := ${OUT_PREFIX}.stl

run-single: ${STL_FILE}

${TOM_FILE}: ${PROG} 
	mkdir -p out
	${PROG} -size ${SIZE} -object ${OBJECT} > ${TOM_FILE}

${STL_FILE}: ${TOM_FILE}
	tomo_to_stl -eps 0.005 < ${TOM_FILE} > ${STL_FILE}
	stmesh_view -format ascii -eps 0.005 -nfGuess ${NF_GUESS} ${STL_FILE} 

endif
# end ${TEST_ARGS} section
######################################################################

clean::

