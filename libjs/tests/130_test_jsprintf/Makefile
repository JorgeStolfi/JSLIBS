# Last edited on 2024-11-18 07:45:03 by stolfi

PROG = test_jsprintf
 
TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

OTHER_C_FLAGS := \
  -Wno-format-zero-length

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test

all: check

check:  do-test

do-test: ${PROG}
	${PROG}

clean:: 
