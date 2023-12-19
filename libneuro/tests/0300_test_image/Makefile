# Last edited on 2021-08-24 16:43:38 by stolfi

PROG := test_neuromat_image

TEST_LIB := libneuro.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

all: check

check:  ${PROG}
	${PROG}
	eom out/*.png

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

clean:: 
	rm -f out/*.{fni,png}
