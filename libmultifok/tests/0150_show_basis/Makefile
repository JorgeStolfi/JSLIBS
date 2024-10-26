# Last edited on 2024-10-11 16:47:26 by stolfi

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
OUT_DIR := out

check:  ${PROG}
	mkdir -p out
	rm -f ${OUT_DIR}/*.ppm ${OUT_DIR}/*.pgm
	./${PROG} \
          -winSize ${NW} \
          -outDir ${OUT_DIR}
	display -title '%f' -filter Box -resize '25600%' ${OUT_DIR}/{weights*-ws,basis*-be}.pgm
