# Last edited on 2023-04-28 20:57:55 by stolfi

PROG = waves_test

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test show
        
FIRST_IMG := out/ppm-400x400/wavys-00.ppm

do-test: ${FIRST_IMG}

${FIRST_IMG}: ${PROG}
	${PROG}

show: ${FIRST_IMG}
	for dir in ppm-400x400 pgm-512x512 pgm-256x256 pgm-128x128; do \
          display -title '%d/%f' -filter Box -resize '200%' `ls out/$${dir}/wavy[us]-*.p?m | sort -r` ; \
        done
        
TEXBANK := ${HOME}/projects/image-collections/texture-bank
        
export: ${FIRST_IMG} 
	for dir in ppm-400x400 pgm-512x512 pgm-256x256 pgm-128x128; do \
          cp -av out/$${dir}/wavy[us]-*.p?m ${TEXBANK}/$${dir}/ ; \
        done

clean::
	/bin/rm -fv out/*/*.{pgm,ppm}
