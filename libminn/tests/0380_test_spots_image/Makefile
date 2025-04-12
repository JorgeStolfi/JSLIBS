# Last edited on 2025-03-31 13:36:55 by stolfi

PROG = test_spots_image 

TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so \
  /usr/lib/x86_64-linux-gnu/libfftw3.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test single-clean single-image all-images

all: 
	./make_all_textures.py
