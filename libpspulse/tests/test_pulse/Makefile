# Last edited on 2020-10-01 20:38:35 by jstolfi

PROG := test_pulse

TEST_LIB := libpspulse.a
TEST_LIB_DIR := ../..

JS_LIBS =  \
  libbezier.a\
  libpspulse.a \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: show

all: check

check:  clean run show

OUT_DIR := out

clean:: 
	-/bin/rm -f ${OUT_DIR}/*.plt ${OUT_DIR}/*.png

run: ${PROG}
	rm -f ${OUT_DIR}/*.plt
	mkdir -p ${OUT_DIR}
	${PROG}

KINDS := B N H

show: ${PROG}
	for kd in ${KINDS}; do \
          for f in `ls ${OUT_DIR}/*-$$kd-*.plt | sort`; do \
            show-derivatives 2 $$f; \
          done; \
        done
  
