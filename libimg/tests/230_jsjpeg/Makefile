# Last edited on 2023-02-27 11:38:23 by stolfi

PROG := test_jsjpeg

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a\
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

VARIANTS := GRY RGB

QUALITIES := 100 099 096 091 084 075 064

QUALVARS := \
  ${addsuffix -GRY,${QUALITIES}} \
  ${addsuffix -RGB,${QUALITIES}}
  
PREFIX := test
 
RDFILES :=  ${addsuffix .jpg,${addprefix out/${PREFIX}-rd-,${VARIANTS}}}
WRFILES :=  ${addsuffix .jpg,${addprefix out/${PREFIX}-wr-,${QUALVARS}}}

all: check

# check:  clean wrfiles
# check:  clean rdfiles
check:  wrfiles

SHOW := display -title '%f' -filter Box -resize '200%'
        
rdfiles: ${RDFILES}

wrfiles: ${PROG} ${RDFILES}
	rm -f ${WRFILES}
	${PROG}
	identify ${RDFILES} ${WRFILES}
	-if [[ -r ${word 1,${WRFILES}} ]]; then ${SHOW} ${WRFILES} ; fi

clean::
	/bin/rm -f ${WRFILES} ${RDFILES}
  
# ====================================================================
# Data files:

# Original images:
ORFILES :=  \
  data/${PREFIX}-orig-GRY.pgm \
  data/${PREFIX}-orig-RGB.ppm

out/%-rd-GRY.jpg: data/%-orig-GRY.pgm
	convert data/$*-orig-GRY.pgm -quality 98 $@
	identify $@
	-if [[ -r $@ ]]; then ${SHOW} $@; fi

out/%-rd-RGB.jpg: data/%-orig-RGB.ppm
	convert data/$*-orig-RGB.ppm -quality 98 $@
	identify $@
	--if [[ -r $@ ]]; then ${SHOW} $@; fi
