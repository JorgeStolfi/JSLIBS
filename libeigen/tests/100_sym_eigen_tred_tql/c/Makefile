# Last edited on 2024-12-05 16:17:30 by stolfi

PROG = test_sym_eigen_tred_tql_C
 
TEST_LIB := libeigen.a
TEST_LIB_DIR := ../../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: show

OUT_DIR := out

all: check

OUT_DIR := out
OUT_FILE_1 := ${OUT_DIR}/version1.txt
OUT_FILE_2 := ${OUT_DIR}/version2.txt

${OUT_FILE_1}: check
${OUT_FILE_2}: check

check: ${PROG} Makefile
	mkdir -p ${OUT_DIR}
	${PROG} > out/.log 2>&1 
	@for ff in ${OUT_FILE_1} ${OUT_FILE_2} ; do \
	  gg="$${ff/out/ref}" ; \
	  echo "=== $${ff} $${gg} ===" ; \
          prdiff -Bb $${ff} $${gg} ; \
        done

show: ${OUT_FILE_1} ${OUT_FILE_2}
	@for ff in ${OUT_FILE_1} ${OUT_FILE_2} ; do \
	  echo "=== $${ff} ===" ; \
	  cat $${ff} ; \
        done

clean::
	-/bin/rm -f ${PROG} *.o ${OUT_FILE_1} ${OUT_FILE_2}
