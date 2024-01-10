# Last edited on 2011-06-06 17:45:33 by stolfi

PROG  := msm_test_220_refining
  
TEST_LIB := libmsmatch.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libimg.a \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run show

all: clean build run show
# all: clean run
# all: clean build

SEQ_LENGTH := 400
DEN = 4

CIRC_X := F
CIRC_Y := F

MUT_DEV  := 0.05
DEL_PROB := 0.10

DELTA :=   10
KAPPA :=   40
MAX_UNP :=  6

# Values are YES or blank:
MAY_EXPAND := YES
MAY_SHRINK := 

BONUS_EQL   := 3.00
PENALTY_DIF := 2.00
PENALTY_SKP := 1.50

N_CANDS := 3
REPEAT := 10

OUTDIR := out

clean::
	-/bin/rm -f ${OUTDIR}/*.pgm ${OUTDIR}/*.ppm
  
run: ${PROG}
	${PROG} \
          -seqLength ${SEQ_LENGTH} -den 4 \
          -seqCirc ${CIRC_X} ${CIRC_Y} \
          -mutDev ${MUT_DEV} -delProb ${DEL_PROB} \
          -delta ${DELTA} -kappa ${KAPPA} -maxUnp ${MAX_UNP} \
          ${MAY_EXPAND:YES=-mayExpand} ${MAY_SHRINK:YES=-mayShrink} \
          -nCands ${N_CANDS} \
          -repeat ${REPEAT} \
          ${OUTDIR}

# PSVIEW := gv
PSVIEW := okular

show: ${OUTDIR}/fin-cd.pgm
	for psf in ${wildcard ${OUTDIR}*.eps} ; do \
          echo "$$psf" ; ${PSVIEW} $$psf ; \
        done
	display \
          -title '%f' \
          -filter Point -resize '<700x700' \
          ${OUTDIR}*-eq.ppm \
          ${OUTDIR}*-ini-cd.pgm \
          ` ( cd ${OUTDIR} && ls [0-9][0-9][0-9][0-9][0-9]-[0-9].ppm ) \
              | sort -t- -k1,1n -k2,2n \
          ` \
          ${OUTDIR}/*fin-cd.pgm
