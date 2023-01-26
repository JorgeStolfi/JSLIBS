# Last edited on 2023-01-22 03:56:32 by stolfi

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

do-test: ${PROG}
	${PROG}

show:
	for dir in ppm-400x400 pgm-512x512 pgm-256x256 pgm-128x128; do \
          display -title '%d/%f' -filter Box -resize '200%' `ls out/$${dir}/wavys-*.p?m | sort` ; \
        done
        
TEXBANK := ${HOME}/projects/image-collections/texture-bank
        
export:    
	for dir in ppm-400x400 pgm-512x512 pgm-256x256 pgm-128x128; do \
          cp -av out/$${dir}/wavys-*.p?m ${TEXBANK}/$${dir}/ ; \
        done

clean::
	/bin/rm -fv ${OUT_PREFIX}-*.{pgm,ppm}
