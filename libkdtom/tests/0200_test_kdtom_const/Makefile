# Last edited on 2021-07-12 05:33:49 by jstolfi

PROG = test_kdtom_const
 
TEST_LIB := libkdtom.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libppv.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

TESTS := \
  .DUMMY

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

# TEST_NAME := tt_${SIZE}_${OBJECT}

# OUT_PREFIX := out/${TEST_NAME}

run-single: ${PROG} 
	mkdir -p out
	${PROG}

endif
# end ${TEST_ARGS} section
######################################################################

clean::

