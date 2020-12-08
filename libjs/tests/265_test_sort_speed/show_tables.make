# Last edited on 2011-05-29 10:48:12 by stolfi 

JOBNAME := show_tables
SIZE    := letter
FMT     := ps

PARTS := \
  test_sort_speed-algs.tex \
  test_sort_speed-ncmp.tex \
  test_sort_speed-time.tex

NOTPARTS :=

MACROS :=

BIBS :=
 
SUBDIRS :=

NEWTEX := export TEXINPUTS=".:${STOLFIHOME}/tex/inputs::"

RUN_LATEX := ${STOLFIHOME}/bin/run-latex
PSVIEW    := okular
PDFVIEW   := acroread

.PHONY: all clean subdirs display-ps display-pdf

all: subdirs ${JOBNAME}.dvi ${JOBNAME}.${FMT} display-${FMT}

subdirs:
	for dir in ${SUBDIRS} ; do \
	  ${MAKE} -f ${JOBNAME}.make SUBDIR=$$dir $$dir.done ; \
	done

######################################################################
# For recursive make - caller must define SUBDIR

SUBDIR := SUBDIR.IS.UNDEFINED
ifneq "${SUBDIR}/" "SUBDIR.IS.UNDEFINED/"

${SUBDIR}.done: ${SUBDIR}
	( cd ${SUBDIR} && ${MAKE} all ) ; /bin/touch ${SUBDIR}.done

endif
# End SUBDIR make
######################################################################

${JOBNAME}.dvi: ${JOBNAME}.tex ${PARTS} ${BIBS} ${SUBDIRS} ${MACROS}
	@echo '= PARTS ='
	@echo ${PARTS} | tr ' ' '\012'
	${RUN_LATEX} ${JOBNAME}.tex
  
${JOBNAME}.ps: ${JOBNAME}.dvi ${SUBDIRS}
	${NEWTEX}; dvips -t ${SIZE} -o ${JOBNAME}.ps -D600 ${JOBNAME}
  
display-ps: ${JOBNAME}.ps
	${PSVIEW} ${JOBNAME}.ps
  
${JOBNAME}.pdf: ${JOBNAME}.ps
	ps-to-pdf < ${JOBNAME}.ps > ${JOBNAME}.pdf
  
display-pdf: ${JOBNAME}.pdf
	${PDFVIEW} ${JOBNAME}.pdf
  
scratch:
	/bin/rm -f ${JOBNAME}.dvi ${JOBNAME}.ps ${JOBNAME}.pdf
  
tidy:
	/bin/rm -f ${JOBNAME}.aux ${JOBNAME}.blg ${JOBNAME}.bbl
	/bin/rm -f ${JOBNAME}.log ${JOBNAME}.errs 
	/bin/rm -f *.dvi~ *.ps~ .pdf~
	/bin/rm -f *.dvi~~ *.ps~~ .pdf~~
	/bin/rm -f *.aux~ *.blg~ *.bbl~ *.log~ *.errs~ core
	/bin/rm -f *.aux~~ *.blg~~ *.bbl~~ *.log~~ *.errs~~
