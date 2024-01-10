# Last edited on 2011-06-06 17:41:36 by stolfi

PROG  := msm_test_120_mapping

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

NCANDS := 5

OUTNAME := out/${CIRC_A}-${CIRC_B}

clean::
	-/bin/rm -f ${OUTNAME}-*.pgm
  
run: ${PROG} 
	${PROG} \
          ${LEN_A} ${CIRC_A} \
          ${LEN_B} ${CIRC_B} \
          ${NCANDS} \
          ${OUTNAME}

show: ${OUTNAME}-00-v.pgm
	display \
          -title '%f' \
          -filter Point -resize '<500x500' \
          ` ls \
              ${OUTNAME}*-[0-9][0-9]-[vm].pgm \
            | sort -t- -k4,4nr -k5,5 \
          `
