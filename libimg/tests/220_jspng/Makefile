# Last edited on 2023-02-26 10:36:30 by stolfi

PROG := test_jspng

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

NOTHING := 
BLANK := ${NOTHING} ${NOTHING}
COMMA := ,

DATADIRS := official stolfi
DATADIRS_COMMAS := ${subst ${BLANK},${COMMA},${strip ${DATADIRS}}}
DATAFILES := ${wildcard 00-DATA/{${DATADIRS_COMMAS}}/*.{png,dir}}

OUTDIRS := ${addprefix out/, ${DATADIRS}}
OUTFILES := ${addsuffix /test-results.dir, ${OUTDIRS}}

all: check

check:  ${PROG} ${DATAFILES}
	@echo "DATADIRS_COMMAS = ${DATADIRS_COMMAS}"
	mkdir -p ${OUTDIRS}
	rm -fv ${OUTFILES}
	${PROG}
	cat ${OUTFILES}

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

clean:: 
	rm -f ${OUTFILES}
