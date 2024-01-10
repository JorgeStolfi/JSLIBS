# Last edited on 2011-06-06 22:32:41 by stolfi

PROG := msm_test_110_cands

TEST_LIB := libmsmatch.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
  
.PHONY:: run show

all: clean run show
# all: clean run

LEN_A := 600
CIRC_A := T

LEN_B := 500
CIRC_B := T

MAXCANDS := 5

OUTNAME := out/${CIRC_A}-${CIRC_B}

OUT_PGM := ${OUTNAME}-cd.pgm
OUT_CDV := ${OUTNAME}-cd.cdv

clean:: 
	-/bin/rm -f ${OUTNAME}-*.pgm
  
run: ${PROG} 
	${PROG} \
          ${LEN_A} ${CIRC_A} \
          ${LEN_B} ${CIRC_B} \
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
