# Last edited on 2019-04-09 22:01:44 by jstolfi

PROG = test_encode_gamma

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

# TEXTURE_OPTS := -dots -vertical
TEXTURE_OPTS :=

OUT_DIR := out

OUTNAME := ${OUT_DIR}/test-enc

# KIND := BT709
KIND := interp
# KIND := generic

ifeq "/${KIND}" "/BT709"
  GAMMA_OPTS := -BT709
endif

ifeq "/${KIND}" "/interp"
  GAMMA_OPTS := -interp \
    0.125 0.330 \
    0.250 0.500 \
    0.500 0.650 \
    0.625 0.690 \
    0.750 0.770 \
    0.875 0.880
endif

ifeq "/${KIND}" "/generic"
  # Gamma and bias values must be in the format "N.NNN"
  # For Sun monitors:
  GAMMA := 0.555
  BIAS  := 0.200

  STEP := 1.200

  GAMMA_OPTS := -generic ${GAMMA} ${BIAS} -step ${STEP}
endif

do-test: ${PROG} ${OUT_DIR}
	${PROG} ${GAMMA_OPTS} ${TEXTURE_OPTS}
	display -title '%f' -gamma 1.000 ${OUTNAME}.pgm
	show_sample_map.sh ${OUTNAME}

clean::
	/bin/rm -fv ${OUT_DIR}/test-*.ppm

${OUT_DIR}:
	mkdir -p ${OUT_DIR}
