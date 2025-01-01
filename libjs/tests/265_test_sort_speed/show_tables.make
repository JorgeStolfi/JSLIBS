# Last edited on 2024-12-31 00:06:10 by stolfi 

JOBNAME := show_tables
SIZE    := letter

PARTS := \
  out/algs.tex \
  out/ncmp.tex \
  out/time.tex

NOTPARTS :=

MACROS :=

BIBS :=
 
SUBDIRS :=

NEWTEX := export TEXINPUTS=".:${STOLFIHOME}/tex/inputs::"

RUN_LATEX := ${STOLFIHOME}/bin/run-pdflatex
PDFVIEW   := evince

.PHONY: all clean-ps display-pdf

PDF_OUT := out/${JOBNAME}.pdf

all: ${PDF_OUT} display-pdf}

${PDF_OUT}: ${JOBNAME}.tex ${PARTS} ${BIBS} ${MACROS}
	@echo 'PARTS = ${PARTS}'
	rm -f ${PDF_OUT}
	cd out && \
          if [[ ! ( -L ${JOBNAME}.tex ) ]]; then \
            ln -s ..$/${JOBNAME}.tex; \
          fi ; \
          {RUN_LATEX} ${JOBNAME}.tex ; \
          rm *.aux *.blg *.bbl *.log *.errs *~
  
display-pdf: ${PDF_OUT)
	${PDFVIEW} ${PDF_OUT}
