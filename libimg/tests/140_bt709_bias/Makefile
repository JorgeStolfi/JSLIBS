# Last edited on 2024-12-20 14:46:58 by stolfi

PROG = test_bt709_bias

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

OPTIONS := -max -weight uniform -expo 0.450 -bias 0.0415

# OPTIONS := -max -weight uniform -expo 0.4501
# OPTIONS := -max -weight uniform -bias 0.0415

do-test: ${PROG}
	${PROG} ${OPTIONS}
	plot_gamma_table.sh 'BT.709' 20 5 'encoding' "${OPTIONS} `cat out/optim.txt`" out/enc
	plot_gamma_table.sh 'BT.709' 20 5 'decoding' "${OPTIONS} `cat out/optim.txt`" out/dec

clean::
	/bin/rm -fv test-*.ppm
