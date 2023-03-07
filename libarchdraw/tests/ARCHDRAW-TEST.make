# Last edited on 2023-02-21 04:39:40 by stolfi
# Makefile for test of {libarchdraw}
# Caller must define ${PROG} ${VERS_TAGS} ${DOC_TAGS} ${EPS_TAGS}

TEST_LIB := libarchdraw.a
TEST_LIB_DIR := ../..

SHELL := /bin/bash

JS_LIBS := \
  libgeo.a \
  libeps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-prog show-doc show-eps show-pts

PREFIX := ${PROG}

IN_PREFIX := in/${PREFIX}
OUT_PREFIX := out/${PREFIX}

OUT_TXT := ${OUT_PREFIX}_P.txt
OUT_TXT_PDF := ${OUT_PREFIX}_P.pdf

all: check

check: run-prog show-doc show-eps show-pts

run-prog: ${OUT_TXT} 

PSVIEW := evince

${OUT_TXT}: ${PROG}
	/bin/rm -f ${OUT_PREFIX}_*.{txt,eps,pdf}
	${PROG}

# Combine page files "${OUT_PREFIX}_${VERSION}_${PAGE}_${DOC}.eps"
# into a single "${OUT_PREFIX}_${VERSION}_${DOC}.pdf":
show-doc: ${OUT_TXT}
	if [[ -s ${OUT_TXT} ]]; then \
          for vers in ${VERS_TAGS} ; do \
            for suff in ${DOC_TAGS} ; do \
              efiles=( `ls ${OUT_PREFIX}_$${vers}_*_$${suff}.eps | sort` ); \
              pfile="${OUT_PREFIX}_$${vers}_$${suff}.pdf" ; \
              merge_eps_to_pdf $${efiles[@]} > $$pfile ; \
              ${PSVIEW} $$pfile ; \
            done ; \
          done ; \
        fi

# Display isolated EPS file "${OUT_PREFIX}_{TAG}.eps":
show-eps: ${OUT_TXT}
	if [[ -s ${OUT_TXT} ]]; then \
          for tag in ${EPS_TAGS} ; do \
            efile=${OUT_PREFIX}_$${tag}.eps ; \
            ${PSVIEW} $$efile ; \
          done ; \
        fi

# Convert list of points to PDF and show:
show-pts: ${OUT_TXT}
	if [[ -s ${OUT_TXT} ]]; then \
	  iso-to-ps-hack -12pt < ${OUT_TXT} | ps2pdf - ${OUT_TXT_PDF} ; \
        fi
	${PSVIEW} ${OUT_TXT_PDF}
