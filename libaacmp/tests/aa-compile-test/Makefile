# Last edited on 2008-01-14 23:29:18 by stolfi

PROG := aa-compile-test
  
JS_LIBS := \
  libaa.a \
  libia.a \
  libflt.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

ENAME := foobar

PCODE := ${ENAME}.pcode
  
EXPR := ${ENAME}.expr
  
CCODE := ${ENAME}.c

IGNORE := ${CCODE}
  
.PHONY:: run
.PRECIOUS: ${CCODE}

all: check

check:  clean run

clean:: 
	-/bin/rm -f ${ENAME}
  
run: ${PROG}
	${PROG}
