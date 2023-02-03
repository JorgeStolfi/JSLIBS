# Last edited on 2023-01-30 21:08:41 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := mf_0180_diff_ops
SCRIPT := run_0180_diff_ops.sh

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a

OTHER_LIBS :=

.PHONY:: check

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 

BASIS_NAME := CANC
TERM_SET := SQR
NOISE := 0.04

all: ${PROG} ${SCRIPT} Makefile
	./${SCRIPT} ${BASIS_NAME} ${TERM_SET} ${NOISE}
	@echo "did it" 1>&2
