# Last edited on 2020-10-11 15:35:57 by jstolfi

PROG := msm_test_110_cands

TEST_LIB := libmsmatch.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libimg.a \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
  
.PHONY:: run show

all: check

check:  clean run show
# check:  clean run

LEN_A := 600

LEN_B := 500

MAXCANDS := 5

OUTDIR := out
OUTNAME := ${OUTDIR}/test

OUT_PGM := ${OUTNAME}-cd.pgm
OUT_CDV := ${OUTNAME}-cd.cdv

clean:: 
	-/bin/rm -f ${OUTDIR}/*.pgm ${OUTDIR}/*.ppm
  
run: ${PROG} 
	mkdir -p ${OUTDIR}
	${PROG} \
          ${LEN_A}  \
          ${LEN_B}  \
          ${MAXCANDS} \
          ${OUTNAME}

show: ${OUT_PGM}
	@echo "=== ${OUT_CDV} =================="
	cat ${OUT_CDV}
	@echo "================================="
	display \
          -title '%f' \
          -filter Point -resize '<500x500' \
          ` ls \
              ${OUT_PGM} \
            | sort -t- -k4,4nr \
          `
