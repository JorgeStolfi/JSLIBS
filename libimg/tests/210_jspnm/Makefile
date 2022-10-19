# Last edited on 2022-10-19 09:12:53 by stolfi

PROG := test_jspnm

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/libnetpbm.so.10

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

VARIANTS := \
  txt-uni.pbm \
  raw-uni.pbm \
  \
  txt-sma.pgm \
  txt-big.pgm \
  raw-sma.pgm \
  raw-big.pgm \
  \
  txt-sma.ppm \
  txt-big.ppm \
  raw-sma.ppm \
  raw-big.ppm
 
PREFIX := test

RDFILES :=  ${addprefix out/${PREFIX}-rd-,${VARIANTS}}
WRFILES :=  ${addprefix out/${PREFIX}-wr-,${VARIANTS}}

all: check

# check:  clean wrfiles
# check:  clean rdfiles
check:  wrfiles

SHOW := display -title '%f%'
        
rdfiles: ${RDFILES}

wrfiles: ${PROG} ${RDFILES}
	${PROG}
	pnmfile ${WRFILES}
	-${SHOW} ${WRFILES}

clean::
	/bin/rm -f ${WRFILES}
  
# ====================================================================
# Data files:

OUTDIR := out

${OUTDIR}/test-rd-raw-uni.pbm: ${OUTDIR} data/test-orig.pgm
	cat data/test-orig.pgm \
	  | pnmdepth 65535 \
	  | pnmscale -pixels 20000 \
	  | pnmnorm -bpercent 0.00001 -wpercent 0.00001 \
	  | pgmtopbm -threshold -value 0.675 \
	  > $@
	pnmfile $@
	-${SHOW} $@

${OUTDIR}/test-rd-txt-uni.%: ${OUTDIR} data/test-orig.pgm
	cat data/test-orig.pgm \
	  | pnmdepth 65535 \
	  | pnmscale -pixels 20000 \
	  | pnmnorm -bpercent 0.00001 -wpercent 0.00001 \
	  | pgmtopbm -threshold -value 0.675 \
	  | pnmnoraw \
	  > $@
	pnmfile $@
	-${SHOW} $@

${OUTDIR}/test-rd-raw-sma.%: ${OUTDIR} data/test-orig.%
	cat data/test-orig.$* \
	  | pnmdepth 65535 \
	  | pnmscale -pixels 20000 \
	  | pnmnorm -bpercent 0.00001 -wpercent 0.00001 \
	  | pnmdepth 255 \
	  > $@
	pnmfile $@
	-${SHOW} $@


${OUTDIR}/test-rd-txt-sma.%: ${OUTDIR} data/test-orig.%
	cat data/test-orig.$* \
	  | pnmdepth 65535 \
	  | pnmscale -pixels 20000 \
	  | pnmnorm -bpercent 0.00001 -wpercent 0.00001 \
	  | pnmdepth 255 \
	  | pnmnoraw \
	  > $@
	pnmfile $@
	-${SHOW} $@

${OUTDIR}/test-rd-raw-big.%: ${OUTDIR} data/test-orig.%
	cat data/test-orig.$* \
	  | pnmdepth 65535 \
	  | pnmscale -pixels 20000 \
	  | pnmnorm -bpercent 0.00001 -wpercent 0.00001 \
	  > $@
	pnmfile $@
	-${SHOW} $@

${OUTDIR}/test-rd-txt-big.%: ${OUTDIR} data/test-orig.%
	cat data/test-orig.$* \
	  | pnmdepth 65535 \
	  | pnmscale -pixels 20000 \
	  | pnmnorm -bpercent 0.00001 -wpercent 0.00001 \
	  | pnmnoraw \
	  > $@
	pnmfile $@
	-${SHOW} $@

${OUTDIR}:
	mkdir -p ${OUTDIR}
