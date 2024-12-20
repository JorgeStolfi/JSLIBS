# Last edited on 2024-12-20 14:18:43 by stolfi

PROG = test_sRGB_bias

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libminu.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

OPTIONS := -max -weight uniform -expo 0.4000 -bias 0.0200

# OPTIONS := -max -weight uniform -expo 0.4000
# OPTIONS := -max -weight uniform -bias 0.0200

do-test: ${PROG}
	${PROG} ${OPTIONS}
	plot_gamma_table.sh 'sRGB' 10 50 'encoding' "${OPTIONS} `cat out/optim.txt`" out/enc
	plot_gamma_table.sh 'sRGB' 10 50 'decoding' "${OPTIONS} `cat out/optim.txt`" out/dec

clean::
	/bin/rm -fv test-*.ppm
