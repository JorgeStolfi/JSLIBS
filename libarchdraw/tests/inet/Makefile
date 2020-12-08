# Last edited on 2012-12-01 15:44:00 by stolfilocal

PROG := inet
TEST_LIB := libarchdraw.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libps.a \
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
PDFVIEW := okular

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
