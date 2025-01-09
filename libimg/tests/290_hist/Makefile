# Last edited on 2024-12-24 10:50:15 by stolfi

PROG = test_float_image_hist

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: clean ${PROG}
	${PROG}
	plot_histograms.sh

clean::
	/bin/rm -fv out/*.{hist,png}
