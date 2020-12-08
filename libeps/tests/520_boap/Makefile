PROG := boap
TEST_LIB := libeps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-prog

OUT_PREFIX := out/porch

OUT_FRONT_EPS := ${OUT_PREFIX}_front_skel.eps

all: check

check:  ${OUT_FRONT_EPS}

${OUT_FRONT_EPS}: run-prog

PSVIEW := evince

run-prog: ${PROG}
	/bin/rm -f ${OUT_PREFIX}_*.{txt,ps,eps,pdf}
	${PROG}
	if [[ -r ${OUT_FRONT_EPS} ]]; then  \
	  ls -l ${OUT_FRONT_EPS}; \
          for efile in ${wildcard ${OUT_PREFIX}*.eps} ; do \
            pfile="$${efile%%.*}.pdf" ; \
            ( echo "$$efile --> $$pfile" 1>&2 ) ; \
            epstopdf $$efile $$pfile ; \
            ${PSVIEW} $$pfile ; \
          done ; \
        fi

# If/when there is a text output:
OUT_TXT := ${OUT_PREFIX}.txt
OUT_TXT_PS := ${OUT_PREFIX}.ps
OUT_TXT_PDF := ${OUT_PREFIX}.pdf

junk:
	iso-to_ps-hack -12pt < ${OUT_TXT} > ${OUT_TXT_PS}
	ps2pdf ${OUT_TXT_PS} ${OUT_TXT_PDF}
	${PSVIEW} ${OUT_TXT_PDF}
