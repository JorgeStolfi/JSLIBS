# Last edited on 2020-10-11 15:28:59 by jstolfi

PROG = test_spectrum_table

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libfftw3.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test do-show do-show-single

TESTS_MISC :=  \
  real:3:0240:0180:0000:0001:0200 \
  peak:1:0007:0007:0001:0002:0100

TESTS_SMALL_WAVE := \
  wave:1:0007:0007:0000:0000:0500 \
  wave:1:0007:0007:0001:0000:0500 \
  wave:1:0007:0007:0001:0001:0500 \
  wave:1:0007:0007:0002:0000:0500 \
  wave:1:0007:0007:0002:0001:0500 \
  wave:1:0007:0007:0002:0002:0500 \
  wave:1:0007:0007:0003:0000:0500 \
  wave:1:0007:0007:0003:0001:0500 \
  wave:1:0007:0007:0003:0002:0500 \
  wave:1:0007:0007:0003:0003:0500
  
TESTS_SMALL_PEAK := \
  peak:1:0005:0005:0001:0002:0005 \
  peak:1:0007:0007:0001:0002:0100 \
  peak:1:0007:0009:0001:0003:0011

TESTS_MEDIUM := \
  bump:3:0072:0045:0024:0022:0045 \
  peak:3:0072:0045:0002:0003:0045 \
  real:3:0232:0188:0000:0000:0045 \
  real:3:0240:0180:0000:0001:0045 \
  real:3:0240:0180:0000:0002:0045 \
  wave:3:0072:0045:0000:0022:0045 \
  wave:3:0072:0045:0024:0015:0045 \
  wave:3:0072:0045:0036:0000:0045 \
  wave:3:0072:0045:0036:0022:0045
  
#TESTS := ${TESTS_SMALL_WAVE}
TESTS := ${TESTS_MISC}

all: check

check:  clean do-test do-show

do-test:  ${PROG}
	for tt in ${TESTS} ; do \
          ${MAKE} TEST_ARGS="$$tt" do-test do-show ; \
        done
        
######################################################################
# For recursive make -- caller must define ${TEST_ARGS}

TEST_ARGS := UNDEFINED
ifneq "/${TEST_ARGS}" "/UNDEFINED"

TEST_FIELDS := ${subst :, ,${TEST_ARGS}}

KIND :=        ${word 1,${TEST_FIELDS}}
CHNS :=        ${word 2,${TEST_FIELDS}}
COLS :=        ${word 3,${TEST_FIELDS}}
ROWS :=        ${word 4,${TEST_FIELDS}}
KX :=          ${word 5,${TEST_FIELDS}}
KY :=          ${word 6,${TEST_FIELDS}}
BINS :=        ${word 7,${TEST_FIELDS}}

# Decide the extension of the image file based on its channel count:
EXT := duh
ifeq "/${CHNS}" "/1"
  EXT := pgm
endif
ifeq "/${CHNS}" "/3"
  EXT := ppm
endif

MAG := 400%

PSVIEW := okular

NAME := out/${KIND}-${COLS}x${ROWS}-${KX}-${KY}

do-test: ${PROG}
	${PROG} ${KIND} ${CHNS} ${COLS} ${ROWS} ${KX} ${KY} ${BINS}

do-show:         
	plot_spectrum_exact "${NAME}-2-tbex"
	@echo ${PSVIEW} ${NAME}-2-tbex-in.eps &
	@echo ${PSVIEW} ${NAME}-2-tbex-ip.eps &
	@echo ${PSVIEW} ${NAME}-2-tbex-cn.eps &
	@echo ${PSVIEW} ${NAME}-2-tbex-cp.eps &
	@echo ${PSVIEW} ${NAME}-2-tbex-cr.eps &
	plot_spectrum_binned "${NAME}-3-tbbn"
	@echo ${PSVIEW} ${NAME}-3-tbbn-b.eps &
	@echo ${PSVIEW} ${NAME}-3-tbbn-n.eps &
	${PSVIEW} ${NAME}-3-tbbn-p.eps &
	@echo ${PSVIEW} ${NAME}-3-tbbn-r.eps &
	plot_spectrum_binned "${NAME}-4-tbcv"
	@echo ${PSVIEW} ${NAME}-4-tbcv-b.eps &
	@echo ${PSVIEW} ${NAME}-4-tbcv-n.eps &
	${PSVIEW} ${NAME}-4-tbcv-p.eps &
	@echo ${PSVIEW} ${NAME}-4-tbcv-r.eps &
	display \
            -title '%f' \
            -filter Point \
            -resize '${MAG}' \
            ${NAME}-0-orig.${EXT} \
            ${NAME}-1-pwrs.${EXT}
	@echo "done."

endif
# End of ${OUT_PPM_FILE} section
######################################################################

clean::
	/bin/rm -fv out/*.{ppm,pgm,txt,eps}
