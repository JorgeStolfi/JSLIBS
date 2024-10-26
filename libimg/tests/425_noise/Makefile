# Last edited on 2024-10-26 05:22:40 by stolfi

PROG = test_noise 

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so \
  /usr/lib/x86_64-linux-gnu/libfftw3.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test show
        
FIRST_IMG := out/png-0640x0480/img-0.000000-0.000000-sqF.png

do-test: ${FIRST_IMG}

${FIRST_IMG}: ${PROG}
	${PROG}

show: ${FIRST_IMG}
	for sr in 0640x0480:200 0400x0400:200 0512x0512:100 0256x0256:200 0128x0128:400 1024x1024:100; do \
          size="$${sr/:*/}" ; mag="$${sr/*:/}"; \
          dir="out/png-$${size}"; \
          display -title "$${size} %f" -filter Box -resize "$${mag}%" \
            `ls $${dir}/{img,pwr}-*.png | sort -t/ -k3.4` ; \
        done
        
TEXBANK := ${HOME}/projects/image-collections/texture-bank
        
export: ${FIRST_IMG} 
	for se in 400:ppm 512:pgm 256:pgm 128:pgm ; do \
	  size="$${se/:*/}" ; ext="$${se/*:/}"; \
          idir="out/png-0$${size}x0$${size}"; \
          odir="${TEXBANK}/$${ext}-$${size}x$${size}"; \
          convert $${odir}/img-0.000000-0.000000-cmF-sqT.png $${odir}/noise-01.$$ext} ; \
          convert $${odir}/img-0.200000-0.200000-cmF-sqT.png $${odir}/noise-02.$$ext} ; \
        done

clean::
	/bin/rm -fv out/*/*.{png,png}
