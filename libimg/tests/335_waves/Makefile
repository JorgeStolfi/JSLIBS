# Last edited on 2025-03-24 15:25:22 by stolfi

PROG = test_waves

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test show
        
FIRST_IMG_NAME := wavys-00
FIRST_IMG := out/pgm-512x512/${FIRST_IMG_NAME}.ppm

do-test: ${FIRST_IMG}

${FIRST_IMG}: ${PROG}
	${PROG}

show: ${FIRST_IMG}
	for dir in ppm-400x400 pgm-512x512 pgm-256x256 pgm-128x128; do \
	  if [[ ( -s $${dir}/${FIRST_IMG_NAME}.ppm ) || ( -s $${dir}/${FIRST_IMG_NAME}.pgm ) ]]; then \
            display -title '%d/%f' -filter Box -resize '200%' `ls out/$${dir}/wavy[us]-*.p?m | sort -r` ; \
          fi \
        done
        
TEXBANK := ${HOME}/projects/image-collections/texture-bank
        
export: ${FIRST_IMG} 
	for dir in ppm-400x400 pgm-512x512 pgm-256x256 pgm-128x128; do \
          cp -av out/$${dir}/wavy[us]-*.p?m ${TEXBANK}/$${dir}/ ; \
        done

clean::
	/bin/rm -fv out/*/*.{pgm,ppm}
