# Last edited on 2023-02-15 12:42:52 by stolfi

PROG = testc
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY::

OUT_DIR := out
OUT_NAME := ${OUT_DIR}/${PROG}

all: check

check:  ${OUT_NAME}.txt

${OUT_NAME}.txt: ${PROG}
	mkdir -p ${OUT_DIR}
	${PROG} > ${OUT_NAME}.txt 2>&1
	@echo "=== ${OUT_NAME}.txt ==="
	cat ${OUT_NAME}.txt

clean::
	-/bin/rm -f ${PROG} *.o ${OUT_NAME}.txt
