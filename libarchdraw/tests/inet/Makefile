# Last edited on 2023-02-04 07:06:38 by stolfi

PROG := inet
TEST_LIB := libarchdraw.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libeps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-prog

OUT_NAME := out/inet
OUT_FILES :=  ${addprefix ${OUT_NAME}, -000001}
PS_FILES :=   ${addsuffix .eps, ${OUT_FILES}}
PDF_FILES :=  ${addsuffix .pdf, ${OUT_FILES}}

all: check

check:  ${PDF_FILES}

${PS_FILES}: run-prog

PSVIEW := gv
PDFVIEW := evince

run-prog: ${PROG}
	${PROG}
	${PSVIEW} ${PS_FILES}

%.pdf: %.ps
	ps2pdf $*.ps $*.pdf
	${PDFVIEW} $*.pdf

%.pdf: %.eps
	ps2pdf $*.eps $*.pdf
	${PDFVIEW} $*.pdf

clean::
	rm -f out/*.ps out/*.eps out/*.pdf out/*.txt
