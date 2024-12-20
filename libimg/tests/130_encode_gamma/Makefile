# Last edited on 2024-12-20 09:59:48 by stolfi

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

all:
	rm -rfv ${OUT_DIR}/test-*
	for kind in BT709 sRGB generic interp ; do \
          ${MAKE} KIND=$${kind} do-test ; \
        done

clean::
	/bin/rm -fv ${OUT_DIR}/test-*.ppm

${OUT_DIR}:
	mkdir -p ${OUT_DIR}
        
MAXVAL := 65535

######################################################################
ifneq "/${KIND}" "/"

ifeq "/${KIND}" "/BT709"
  GAMMA_OPTS := -kind BT709
endif

ifeq "/${KIND}" "/sRGB"
  GAMMA_OPTS := -kind sRGB
endif

ifeq "/${KIND}" "/interp"
  GAMMA_OPTS := -kind interp \
    -points \
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
  EXPO := 0.555
  BIAS := 0.200
  GAMMA_OPTS := -kind generic -expo ${EXPO} -bias ${BIAS}
endif

OUTNAME := test-${KIND}
OUTFILE := ${OUT_DIR}/${OUTNAME}.txt

do-test: ${PROG} ${OUT_DIR}
	${PROG} -maxval ${MAXVAL} ${GAMMA_OPTS} 
	show_sample_map.sh ${OUTFILE} ${KIND}

endif
######################################################################
