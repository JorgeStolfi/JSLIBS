# Last edited on 2023-02-04 08:07:31 by stolfi

PROG := test_calib_opt
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

.PHONY:: do-all-tests do-test

JS_LIBS := \
  libgeo.a \
  libeps.a \
  libjs.a
  
TNAMES := ${shell cd in && ls -d [a-z][0-9][0-9]}

all: check

check:  build do-all-tests

do-all-tests: 
	for tn in ${TNAMES} ; do \
	  ${MAKE} TNAME=$${tn} do-test ; \
        done

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

include ../../../../GENERIC-LIB-TEST.make

######################################################################
# Recursive make - client must define ${TNAME}
ifneq "/${TNAME}" "/"

do-test: ${PROG} in/${TNAME} run_test.sh
	./run_test.sh ${TNAME}

endif
# end ${TNAME} section
######################################################################
