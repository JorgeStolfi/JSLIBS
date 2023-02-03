# Last edited on 2023-01-30 12:17:55 by stolfi

PROG := ../../make_CANC_HART_terms.gawk

BASIS_NAME := CANC

ifeq "/${BASIS_NAME}" "/CANC" 
  BEL_NAMES := Soo Smo Spo Som Sop Smm Spp Smp Spm
else ifeq "/${BASIS_NAME}" "/HART" 
  BEL_NAMES := Hoo Hmo Hpo Hom Hop Hmm Hpp Hmp Hpm
endif

OUT_PREFIX := out/${BASIS_NAME}

all: ${PROG} Makefile
	${PROG} \
	  -v belNames="${BEL_NAMES}" \
	  -v outPrefix="${OUT_PREFIX}" 
	@echo "=== terms ===" 1>&2
	@cat ${OUT_PREFIX}-terms.txt 1>&2
	@echo "=== tix ===" 1>&2
	@cat ${OUT_PREFIX}-tix.txt 1>&2
