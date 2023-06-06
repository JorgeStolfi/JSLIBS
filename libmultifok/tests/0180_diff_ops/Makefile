# Last edited on 2023-04-29 04:48:41 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_diff_ops
SCRIPT := run_mfok_diff_ops.sh

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a

OTHER_LIBS :=

.PHONY:: check

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

IMAGE_SET := SMUF
BASIS_NAME := CANC
TERM_SET := ALL
NOISE := 0.02

all: ${PROG} ${SCRIPT} Makefile
	./${SCRIPT} ${IMAGE_SET} ${BASIS_NAME} ${TERM_SET} ${NOISE}
	@echo "did it" 1>&2
