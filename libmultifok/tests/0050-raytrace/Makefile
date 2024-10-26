# Last edited on 2024-10-21 20:20:46 by stolfi

SHELL := /bin/bash

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_raytrace

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: \
  run-prog

all: run-prog show-rays

# Image dimensions:
# IMAGE_NX := 640
# IMAGE_NY := 480

# IMAGE_NX := 160
# IMAGE_NY := 120

IMAGE_NX := 320
IMAGE_NY := 240

IMAGE_TYPE := bullsqr
PIX_SAMPLING := 3
DIR_SAMPLING := 10

SIZE_TAG := ${shell printf "%04dx%04d" ${IMAGE_NX} ${IMAGE_NY}}
SAMPLING_TAG := ${shell printf "hs%02d-kr%02d" ${PIX_SAMPLING} ${DIR_SAMPLING}}
RUN_DIR := out/img-${IMAGE_TYPE}-${SIZE_TAG}-${SAMPLING_TAG}

run-prog: ${PROG} ${RUN_SCRIPT} Makefile
	rm -f ${RUN_DIR}-zf*.png
	./${PROG} \
          -imageSize ${IMAGE_NX} ${IMAGE_NY} \
          -imageType ${IMAGE_TYPE} \
          -pixSampling ${PIX_SAMPLING} \
          -dirSampling ${DIR_SAMPLING}
	display -title '%f' ${RUN_DIR}/zf*/img-blur.png 

show-rays: mfok_plot_rays.sh ${RUN_DIR}
	for frameDir in `cd ${RUN_DIR} && ls -d zf*-df* | sort` ; do \
	  mfok_plot_rays.sh ${RUN_DIR}/$${frameDir}/pixel-rays-0010-0017.txt ; \
        done
