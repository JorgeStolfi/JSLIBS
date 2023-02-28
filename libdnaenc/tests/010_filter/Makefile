# Last edited on 2023-02-27 19:46:36 by stolfi

TEST_LIB := libdnaenc.a
TEST_LIB_DIR := ../..
PROG := test_dnae_filter

JS_LIBS := \
  libmsmatch.a \
  libgeo.a \
  libeps.a \
  libimg.a \
  libjs.a

OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libfftw3.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

DATA_DIR := ${STOLFIHOME}/programs/c/DNA/DATA

SOME_SEQS := \
  small-old:SM_001542 \
  small/bas:250B \
  small/bas:1500A \
  small/bas:2500A

SMALL_SEQS := \
  small/bas:100A \
  small/bas:250A \
  small/bas:250B \
  small/bas:250E \
  small/bas:250E-in01 \
  small/bas:250E-in02 \
  small/bas:250E-in04 \
  small/bas:250E-in08 \
  small/bas:250E-mt01 \
  small/bas:250E-mt02 \
  small/bas:250E-mt04 \
  small/bas:250E-mt08 \
  small/bas:500A \
  small/bas:500A-in01 \
  small/bas:500A-in02 \
  small/bas:500A-in04 \
  small/bas:500A-in08 \
  small/bas:500A-mt01 \
  small/bas:500A-mt02 \
  small/bas:500A-mt04 \
  small/bas:500A-mt08 \
  small/bas:500B \
  small/bas:1500A \
  small/bas:2500A \
  small/bas:2500B \
  small/bas:2500C

# SEQS := small-old/bas:SM_001542
# SEQS := small-old/bas:QQ_002489_s1

# SEQS := ${SMALL_SEQS}
SEQS := ${SOME_SEQS}

MAX_LEVEL := 11

PLOT_WIDTH := 72.0
PLOT_HEIGHT := 18.0
FONT_SIZE := 10.0

INIT_WEIGHTS := wg020
INCR_WEIGHTS := wg060

SHOW := 0
VIEW := 1

# ACTIONS := clean-single run-single show-single view-single
ACTIONS := clean-single run-single

all: check

check:  clean run-all

run-all: 
	for seq in ${SEQS} ; do \
            ${MAKE} SEQ=$$seq clean-single ; \
            ${MAKE} SEQ=$$seq run-single ; \
            if [[ ${SHOW} -ne 0 ]]; then ${MAKE} SEQ=$$seq show-single ; fi ; \
            if [[ ${VIEW} -ne 0 ]]; then ${MAKE} SEQ=$$seq view-single ; fi ; \
        done
     
# ======================================================================
# Recursive make - must define ${SEQ}
ifneq "/${SEQ}" "/"        

ARGS := ${subst :, ,${SEQ}}
IN_DIR :=  ${word 1, ${ARGS}}
IN_NAME := ${word 2, ${ARGS}}

IN_PREFIX := ${DATA_DIR}/${IN_DIR}/${IN_NAME}

OUT_PREFIX := out/${IN_NAME}-${INIT_WEIGHTS}-${INCR_WEIGHTS}

SHOWPS := gv
VIEWEQS := dna_seq_view -resample 0.25 -seqFile
  
run-single: ${PROG}
	${PROG} \
          -maxLevel ${MAX_LEVEL} \
          -initFilter ${shell cat ${DATA_DIR}/weights/${INIT_WEIGHTS}.txt | sed -e 's:[\#].*$$::g'} \
          -initStep 1 \
          -incrFilter ${shell cat ${DATA_DIR}/weights/${INCR_WEIGHTS}.txt | sed -e 's:[\#].*$$::g'} \
          -plotSignals \
          -plotSpectra \
          -plotSize ${PLOT_WIDTH} ${PLOT_HEIGHT} \
          -fontSize ${FONT_SIZE} \
          ${IN_PREFIX} \
          ${OUT_PREFIX}

EPS_FILE_00 := ${OUT_PREFIX}-00.eps
EPS_FILES := ${sort ${wildcard ${OUT_PREFIX}-??.eps} ${wildcard ${OUT_PREFIX}-??-f?.eps}}

EQS_FILE_00 := ${OUT_PREFIX}-00.eqs
EQS_FILES := ${sort ${wildcard ${OUT_PREFIX}-??.eqs}}

show-single: ${EPS_FILE_00}
	@for psf in ${EPS_FILES} ; do \
          echo "$$psf" ; ${SHOWPS} -title $${psf} $${psf} ;  \
        done

view-single: ${EQS_FILE_00}
	@for eqf in ${EQS_FILES} ; do \
          echo "$$eqf" ; ${VIEWEQS} $$eqf ;  \
        done
  
clean-single:: 
	-/bin/rm -f ${OUT_PREFIX}-*.pgm ;
	-/bin/rm -f ${OUT_PREFIX}-*.egs ;
	-/bin/rm -f ${OUT_PREFIX}-*.eps ;
 
endif
# ======================================================================
