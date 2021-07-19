# Last edited on 2021-07-19 04:12:06 by jstolfi

PROG := test_calib_opt
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS := \
  liblmdif.a \
  libgeo.a \
  libps.a \
  libjs.a

all: check

check:  build run-tests

run-all-tests: \
  ${MAKE} TNAME=t01 run-test

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

include ../../../../GENERIC-LIB-TEST.make

######################################################################
# Recursive make - client must define ${TNAME}
ifneq "/${TNAME}" "/"

run-test: ${PROG} in/${TNAME}
	./run_test.sh ${TNAME}

endif
# end ${TNAME} section
######################################################################
