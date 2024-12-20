# Last edited on 2024-12-20 10:39:50 by stolfi

PROG = test_gamma_image

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

all:
	rm -rf ${OUT_DIR}/test-*
	for kind in BT709 sRGB generic interp ; do \
          ${MAKE} KIND=$${kind} do-test ; \
        done

clean::
	/bin/rm -fv ${OUT_DIR}/test-*.ppm

${OUT_DIR}:
	mkdir -p ${OUT_DIR}

######################################################################
ifneq "/${KIND}" "/"

ifeq "/${KIND}" "/generic"
  # Gamma and bias values must be in the format "N.NNN"
  # For Sun monitors:
  EXPO := 0.555
  BIAS  := 0.200
  STEP := 1.200
  GENERIC_OPTS := -expo ${EXPO} -bias ${BIAS} -step ${STEP}
endif

OUTNAME := test-${KIND}
OUTPGM := ${OUT_DIR}/${OUTNAME}.pgm

do-test: ${PROG} ${OUT_DIR}
	${PROG} -kind ${KIND} ${GENERIC_OPTS} ${TEXTURE_OPTS}
	display -title '%f' -gamma 1.000 ${OUTPGM}

endif
######################################################################
