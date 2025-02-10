# Last edited on 2025-01-30 10:38:02 by stolfi

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
OUT_FOLDER := out

check:  ${PROG}
	mkdir -p out
	rm -f ${OUT_FOLDER}/*.ppm ${OUT_FOLDER}/*.pgm
	./${PROG} \
          -winSize ${NW} \
          -outFolder ${OUT_FOLDER}
	display -title '%f' -filter Box -resize '25600%' ${OUT_FOLDER}/{weights*-ws,basis*-be}.pgm
