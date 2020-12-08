# Last edited on 2017-06-15 16:43:10 by stolfilocal

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

# OPTIONS := -mean -weight uniform -gamma 0.450

# OPTIONS := -max -weight uniform -gamma 0.4545454545454545

# OPTIONS := -max -weight uniform -bias 0.0333333333333333
# OPTIONS := -max -weight uniform -bias 0.0390

OPTIONS := -mean -weight uniform -gamma 0.450 -bias 0.0327

do-test: ${PROG}
	${PROG} ${OPTIONS}
	plot_gamma_table.sh 'encoding' "${OPTIONS} `cat out/optim.txt`" out/enc
	plot_gamma_table.sh 'decoding' "${OPTIONS} `cat out/optim.txt`" out/dec

clean::
	/bin/rm -fv test-*.ppm
