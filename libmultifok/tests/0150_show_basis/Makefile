# Last edited on 2023-04-28 17:21:35 by stolfi

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_show_basis

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

NW := 3
OUT_PREFIX := out/bas

check:  ${PROG}
	mkdir -p out
	rm -f ${OUT_PREFIX}*.ppm ${OUT_PREFIX}*.pgm
	./${PROG} \
          -winSize ${NW} \
          -outPrefix ${OUT_PREFIX}
	display -title '%f' -filter Box -resize '12800%' ${OUT_PREFIX}*-{ws,be}.pgm
