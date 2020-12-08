PROG := boas
TEST_LIB := libarchdraw.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libeps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-prog

PREFIX := ${PROG}

IN_PREFIX := in/${PREFIX}
OUT_PREFIX := out/${PREFIX}

OUT_TXT := ${OUT_PREFIX}_p.txt
OUT_PS := ${OUT_PREFIX}_p.ps
OUT_PDF := ${OUT_PREFIX}_p.pdf

all: check

check:  ${OUT_TXT}

${OUT_TXT}: run-prog

PSVIEW := evince

run-prog: ${PROG}
	/bin/rm -f ${OUT_PREFIX}_*.{txt,ps,eps,pdf}
	${PROG}
	if [[ -r ${OUT_TXT} ]]; then  \
          for efile in ${wildcard ${OUT_PREFIX}_*.eps} ; do \
            pfile="$${efile##.*}.pdf" ; \
            ( echo "$$efile --> $$pfile" 1>&2 ) ; \
            epstopdf $$efile $$pfile ; \
            ${PSVIEW} $$pfile ; \
          done ; \
        fi
	iso-to_ps-hack -12pt < ${OUT_TXT} > ${OUT_PS}
	ps2pdf ${OUT_PS} ${OUT_PDF}
	${PSVIEW} ${OUT_PDF}
