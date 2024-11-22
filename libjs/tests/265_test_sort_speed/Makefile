# Last edited on 2024-11-18 06:08:10 by stolfi

PROG := test_sort_speed

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: clean check

check:  do-test

do-test: ${PROG}
	time ${PROG}

clean::
	-/bin/rm -f out/*.tex
