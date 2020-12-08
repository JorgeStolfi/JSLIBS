PROG := miguec
TEST_LIB := libeps.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-prog

OUT_PREFIX := out/control_diagram

OUT_DIAG_EPS := ${OUT_PREFIX}.eps
OUT_DIAG_PDF := ${OUT_PREFIX}.pdf

all: check

check:  ${OUT_DIAG_EPS}

${OUT_DIAG_EPS}: run-prog

PSVIEW := evince

run-prog: ${PROG}
	/bin/rm -f ${OUT_PREFIX}_*.{txt,ps,eps,pdf}
	${PROG}
	if [[ -r ${OUT_DIAG_EPS} ]]; then  \
          for efile in ${OUT_DIAG_EPS} ; do \
            pfile="$${efile%%.*}.pdf" ; \
            ( echo "$$efile --> $$pfile" 1>&2 ) ; \
            epstopdf $$efile $$pfile ; \
            ${PSVIEW} $$pfile ; \
          done ; \
        fi
