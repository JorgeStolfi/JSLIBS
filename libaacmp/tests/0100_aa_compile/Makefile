# Last edited on 2024-11-08 11:16:05 by stolfi

PROG := test_aa_compile

TEST_LIB_DIR := ../..
TEST_LIB := libaacmp.a 
  
JS_LIBS := \
  libaa.a \
  libia.a \
  libflt.a \
  libjs.a

all: check

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

ENAME := foobar
  
EXPR := in/${ENAME}.expr

PCODE := out/${ENAME}.pcode
  
CCODE := out/${ENAME}.c
  
.PHONY:: run

check:  clean run

clean:: 
	-/bin/rm -f out/*.pcode out/*.c
  
run: ${CCODE}

${CCODE}: ${PCODE} ${PROG} Makefile
	${PROG} ${ENAME}
	@echo "......................................................................" 1>&2
	@cat ${CCODE} 1>&2
	@echo "......................................................................" 1>&2
        
ETOP := ${STOLFIHOME}/bin/expr_to_pcode.sh

${PCODE}: ${ETOP} ${EXPR} Makefile
	${ETOP} < ${EXPR} > ${PCODE}
	@echo "......................................................................" 1>&2
	@cat ${PCODE} 1>&2
	@echo "......................................................................" 1>&2
