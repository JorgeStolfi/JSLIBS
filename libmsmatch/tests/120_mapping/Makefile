# Last edited on 2020-10-11 15:32:21 by jstolfi

PROG  := msm_test_120_mapping

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

NCANDS := 5

OUTDIR := out
OUTNAME := ${OUTDIR}/test

clean::
	-/bin/rm -f ${OUTDIR}/*.pgm ${OUTDIR}/*.ppm
  
run: ${PROG} ${OUTDIR}
	mkdir -p ${OUTDIR}
	${PROG} \
          ${LEN_A}  \
          ${LEN_B}  \
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

${OUTDIR}:
	mkdir -p ${OUTDIR}
          
