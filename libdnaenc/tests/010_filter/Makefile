# Last edited on 2020-10-11 15:40:39 by jstolfi

TEST_LIB := libdnaenc.a
TEST_LIB_DIR := ../..
PROG := test_dnae_filter

JS_LIBS := \
  libmsmatch.a \
  libgeo.a \
  libps.a \
  libimg.a \
  libjs.a

OTHER_LIBS := \
  /lib64/libfftw3.so.3

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

DATA_DIR := ${STOLFIHOME}/programs/c/DNA/DATA
IN_DIR := ${DATA_DIR}/small/bas

SMALL_SEQS := \
  100A \
  250A \
  250B \
  250E \
  250E-in01 \
  250E-in02 \
  250E-in04 \
  250E-in08 \
  250E-mt01 \
  250E-mt02 \
  250E-mt04 \
  250E-mt08 \
  500A \
  500A-in01 \
  500A-in02 \
  500A-in04 \
  500A-in08 \
  500A-mt01 \
  500A-mt02 \
  500A-mt04 \
  500A-mt08 \
  500B \
  1500A \
  2500A \
  2500B \
  2500C

# SEQS := SM_001542
# SEQS := QQ_002489_s1
# SEQS := /Dsimulans/Dsimulans_4-norm
# SEQS := AF_018985-droso-1-norm

SEQS := ${SMALL_SEQS}

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

IN_NAME := ${IN_DIR}/${SEQ}

OUT_NAME := out/${SEQ}-${INIT_WEIGHTS}-${INCR_WEIGHTS}

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
          ${IN_NAME} \
          ${OUT_NAME}

EPS_FILE_00 := ${OUT_NAME}-00.eps
EPS_FILES := ${sort ${wildcard ${OUT_NAME}-??.eps} ${wildcard ${OUT_NAME}-??-f?.eps}}

EQS_FILE_00 := ${OUT_NAME}-00.eqs
EQS_FILES := ${sort ${wildcard ${OUT_NAME}-??.eqs}}

show-single: ${EPS_FILE_00}
	@for psf in ${EPS_FILES} ; do \
          echo "$$psf" ; ${SHOWPS} -title $${psf} $${psf} ;  \
        done

view-single: ${EQS_FILE_00}
	@for eqf in ${EQS_FILES} ; do \
          echo "$$eqf" ; ${VIEWEQS} $$eqf ;  \
        done
  
clean-single:: 
	-/bin/rm -f ${OUT_NAME}-*.pgm ;
	-/bin/rm -f ${OUT_NAME}-*.egs ;
	-/bin/rm -f ${OUT_NAME}-*.eps ;
 
endif
# ======================================================================
