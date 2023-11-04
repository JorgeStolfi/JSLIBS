# Last edited on 2023-10-22 03:29:15 by stolfi

TEST_LIB := libneuro.a
TEST_LIB_DIR := ../..
PROG := test_butterworth

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libfftw3.so \

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  ${PROG}
	mkdir -p out
	./${PROG} > out/butfil.txt
	./plot_filter.sh out/butfil.txt
