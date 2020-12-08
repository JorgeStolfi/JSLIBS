# Last edited on 2008-01-15 05:41:35 by stolfi

PROG := aavalid
  
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
