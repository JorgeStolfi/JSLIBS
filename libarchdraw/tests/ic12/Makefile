# Last edited on 2020-10-11 15:22:41 by jstolfi

PROG := ic12
TEST_LIB := libarchdraw.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libps.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-prog export

PREFIX := ${PROG}

IN_PREFIX := in/${PREFIX}
OUT_PREFIX := out/${PREFIX}

OUT_FILES :=  ${addprefix ${OUT_PREFIX}, -d -p}
PS_FILES :=   ${addsuffix .ps, ${OUT_FILES}}
PDF_FILES :=  ${addsuffix .pdf, ${OUT_FILES}}

all: check

check:  ${PDF_FILES}
# check:  ${PS_FILES} export

${PS_FILES}: run-prog

PSVIEW := evince

run-prog: ${PROG}
	/bin/rm -f ${OUT_PREFIX}-d.ps ${wildcard ${OUT_PREFIX}-*.eps}
	${PROG}
	if [[ -r ${OUT_PREFIX}-d.ps ]]; then  ${PSVIEW} ${OUT_PREFIX}-d.ps ; fi
	# for ff in ${OUT_PREFIX}-*.eps ; do ${PSVIEW} $$ff; done
	-display ${OUT_PREFIX}-*.eps
	iso-to-ps-hack -12pt < ${OUT_PREFIX}-p.txt > ${OUT_PREFIX}-p.ps

%.pdf: %.ps
	ps2pdf $*.ps $*.pdf

EXPORT_DIR := ${STOLFIHOME}/bureau/predios/ic-4a/figures

export: run-prog
	if [[ -r ${OUT_PREFIX}-000001.eps ]]; then cp -av ${OUT_PREFIX}-*.eps ${EXPORT_DIR} ; fi
