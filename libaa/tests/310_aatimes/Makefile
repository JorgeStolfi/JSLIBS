# Last edited on 2008-01-14 22:43:44 by stolfi

PROG := aatimes
  
TEST_LIB := libaa.a 
TEST_LIB_DIR := ../..

JS_LIBS := \
  libia.a \
  libflt.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  run

run: ${PROG}
	${PROG}
