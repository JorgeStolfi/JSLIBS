# Last edited on 2023-06-24 20:12:20 by stolfi

PROG := test_drtree_plot

TEST_LIB := libdrtree.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libeps.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test show-epsf

all: check

check:  do-test

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/tdp

MAIN_EPS := ${OUT_PREFIX}-plot.eps

PSVIEW := evince
        
do-test: ${MAIN_EPS} show-epsf
  
${MAIN_EPS} : ${PROG} ${OUT_DIR}
	./${PROG} \
          -nIndivs 300 \
          -nRoots 10 \
          -tStart 100 \
          -tStop 350 \
          -orphans T \
          -ageMax 25 \
          -nchMax 6 \
          -tRef 200 \
          -outPrefix  ${OUT_PREFIX}

show-epsf: ${MAIN_EPS}
	for f in `ls ${OUT_PREFIX}-*.eps | sort` ; do \
	  ${PSVIEW} $$f ; \
        done
        
${OUT_DIR}:
	mkdir -p ${OUT_DIR}

clean::
	-/bin/rm -f ${OUT_PREFIX}-*.eps
  
