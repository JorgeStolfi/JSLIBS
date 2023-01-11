# Last edited on 2023-01-10 20:19:46 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := mf_0100_make_image

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

NX := 160
NY := 120
NR := 400
FOC_DEPTH := 2.0
DENSE := F

OUT_PREFIX := out/test-d${DENSE}

check:  ${PROG}
	mkdir -p out
	rm -f ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.ppm
	./${PROG} \
          -imgSize ${NX} ${NY} \
          -rays ${NR} \
          -dense ${DENSE} \
          -focDepth ${FOC_DEPTH} \
          -zRange 0.0 10.0 -zStep 1.0 \
          -rPix 0.50 \
          -outPrefix ${OUT_PREFIX}
	display -filter Box -resize '400%' ${OUT_PREFIX}*.pgm ${OUT_PREFIX}*.ppm
