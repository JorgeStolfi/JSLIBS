# Last edited on 2013-06-14 01:57:23 by stolfilocal

TEST_LIB := libneuro.a
TEST_LIB_DIR := ../..
PROG := test_butterworth

JS_LIBS := \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib64/libfftw3.so \

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check:  ${PROG}
	mkdir -p out
	./${PROG} > out/butfil.txt
	./plot_filter.sh out/butfil.txt
