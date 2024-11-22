# Last edited on 2024-11-21 21:22:35 by stolfi

PROG = test_sym_eigen
 
TEST_LIB := libgeo.a
TEST_LIB_DIR := ../../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: show

OUT_DIR := out
OUT_NAME := ${OUT_DIR}/${PROG}

all: check show

check:  ${OUT_NAME}.txt

# DAMN MINUS ZEROS!!!!

${OUT_NAME}.txt: ${PROG} Makefile
	mkdir -p ${OUT_DIR}
	( ${PROG} 2>&1 ) | sed -e 's:[-]\(0[.]000*\b\): \1:g' > ${OUT_NAME}.txt

show: ${OUT_NAME}.txt
	@echo "=== ${OUT_NAME}.txt ==="
	cat ${OUT_NAME}.txt

clean::
	-/bin/rm -f ${PROG} *.o ${OUT_NAME}.txt
