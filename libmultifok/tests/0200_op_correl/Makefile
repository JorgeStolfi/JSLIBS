# Last edited on 2023-01-09 12:56:27 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := mf_0200_op_correl

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

NW := 3
NOISE := 0.002
Z_FOCUS := 005.0000
IN_PREFIX := in/test-dF-fd002.0000
OUT_PREFIX := out/test
BASIS_TYPE := D
SQR_ONLY := 1

check:  ${PROG}
	mkdir -p out
	rm -f out/*.ppm out/*.pgm
	./${PROG} \
            -inPrefix ${IN_PREFIX} \
            -winSize ${NW} \
            -basisType ${BASIS_TYPE} \
            -zFocus ${Z_FOCUS} \
            -noise ${NOISE} \
            -outPrefix ${OUT_PREFIX}
	display -filter Box -resize '400%' ${OUT_PREFIX}*.pgm
	do_regression.sh ${OUT_PREFIX}-z${Z_FOCUS} ${BASIS_TYPE} ${NW} ${SQR_ONLY}
	plot_regression.sh ${OUT_PREFIX}-z${Z_FOCUS} ${BASIS_TYPE} ${NW}
            
