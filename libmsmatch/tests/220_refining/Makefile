# Last edited on 2023-02-04 07:05:47 by stolfi

PROG  := msm_test_220_refining
  
TEST_LIB := libmsmatch.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libimg.a \
  libeps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run show

all: check

check:  clean build run show
# check:  clean run
# check:  clean build

SEQ_LENGTH := 400
EXP_SUB = 2

MUT_DEV  := 0.05
DEL_PROB := 0.10

DELTA :=   10
KAPPA :=   40
MAX_UNP :=  6

# Candidate length variation:
EXPAND := 15
SHRINK := 15

BONUS_EQL   := 3.00
PENALTY_DIF := 2.00
PENALTY_SKP := 1.50

N_CANDS := 3
REPEAT := 10

OUTDIR := out
OUTNAME := ${OUTDIR}/test

clean::
	-/bin/rm -f ${OUTDIR}/*.pgm ${OUTDIR}/*.ppm
  
run: ${PROG}
	mkdir -p ${OUTDIR}
	${PROG} \
          -seqLength ${SEQ_LENGTH} \
          -expSub ${EXP_SUB} \
          -mutDev ${MUT_DEV} -delProb ${DEL_PROB} \
          -delta ${DELTA} -kappa ${KAPPA} -maxUnp ${MAX_UNP} \
          -expand ${EXPAND} -shrink ${SHRINK} \
          -nCands ${N_CANDS} \
          -repeat ${REPEAT} \
          ${OUTNAME}

# PSVIEW := gv
PSVIEW := evince

show: ${OUTNAME}-fin-cd.pgm
	for psf in ${wildcard ${OUTDIR}*.eps} ; do \
          echo "$$psf" ; ${PSVIEW} $$psf ; \
        done
	display \
          -title '%f' \
          -filter Point -resize '<700x700' \
          -filter Lanczos -resize '>900x900' \
          ${OUTNAME}*-eq.ppm \
          ${OUTNAME}*-ini-cd.pgm \
          ` ( ls ${OUTNAME}-[0-9][0-9][0-9][0-9][0-9]-[0-9].ppm ) \
              | sort -t- -k1,1n -k2,2n \
          ` \
          ${OUTNAME}*-fin-cd.pgm
