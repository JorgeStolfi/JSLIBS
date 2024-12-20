# Last edited on 2024-12-04 11:28:47 by stolfi

PROG = test_sym_eigen_tridiag
 
TEST_LIB := libeigen.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: show

OUT_DIR := out
OUT_NAME := ${OUT_DIR}/testc

all: check show

OUT_DIR := out
OUT_FILE_1 := ${OUT_DIR}/test1.txt
OUT_FILE_2 := ${OUT_DIR}/test2.txt
OUT_FILE_3 := ${OUT_DIR}/test3.txt

${OUT_FILE_1}: check
${OUT_FILE_2}: check
${OUT_FILE_3}: check

check: ${PROG} Makefile
	mkdir -p ${OUT_DIR}
	${PROG} > out/.log 2>&1 
	@for ff in ${OUT_FILE_1} ${OUT_FILE_2} ${OUT_FILE_3} ; do \
	  gg="$${ff/out/ref}" ; \
	  echo "######### $${ff} $${gg} #########" ; \
          prdiff -Bb $${ff} $${gg} ; \
        done

show: ${OUT_FILE_0} ${OUT_FILE_1} ${OUT_FILE_2}
	@for ff in ${OUT_FILE_1} ${OUT_FILE_2} ${OUT_FILE_3} ; do \
	  echo "######### $${ff} #########" ; \
	  cat $${ff} ; \
        done

clean::
	-/bin/rm -f ${PROG} *.o ${OUT_FILE_1} ${OUT_FILE_2} ${OUT_FILE_3}
