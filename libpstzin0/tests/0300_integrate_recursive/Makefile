# Last edited on 2025-03-16 12:09:52 by stolfi

PROG = test_integrate_recursive

TEST_LIB := libpstzin0.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a \
  libjsextra.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-tests

all: check

# check: do-tests

check: do-gustavo

# FUNC := 01-ramp10-N 
# FUNC := 02-ramp01-N 
# FUNC := 03-ramp11-N 
# FUNC := 04-parabo-N
# FUNC := 20-cpiece-N
# FUNC := 24-cplatv-N    
FUNC := 18-cbramp-N
# FUNC := 21-cplat3-N

# SIZE := 0001x0001
# SIZE := 0002x0002
# SIZE := 0016x0012
SIZE := 0032x0024
# SIZE := 0256x0192

HINTS := 0
REFERENCE := 1
INITIAL_OPTION := zero
INITIAL_NOISE := 0.25

# MAX_LEVEL := 999
MAX_LEVEL := 0

do-tests: ${PROG}
	mkdir -p out
	./run_test.sh ${FUNC} ${SIZE} ${HINTS} ${REFERENCE} ${INITIAL_OPTION} ${INITIAL_NOISE} ${MAX_LEVEL}
 
GU_SCENE := manyob
GU_SIZE := 0512x0422


GU_INPUT_DTSEQ := 2025-03-03-00
GU_INPUT_METHOD := WA


GU_HINTS_DTSEQ := NONE
GU_HINTS_METHOD := MF

# GU_HINTS_DTSEQ := 2025-03-05-00
# GU_HINTS_METHOD := MF

# GU_HINTS_DTSEQ := 2025-03-03-00
# GU_HINTS_METHOD := SY


GU_REFERENCE_DTSEQ := 2025-03-03-00
GU_REFERENCE_METHOD := SY

GU_INITIAL_OPT := zero
GU_INITIAL_NOISE := 0.0

GU_MAXLEVEL := 4

do-gustavo:
	mkdir -p out
	./run_gustavo.sh \
          ${GU_SCENE} ${GU_SIZE} \
          ${GU_INPUT_DTSEQ} ${GU_INPUT_METHOD} \
          ${GU_HINTS_DTSEQ} ${GU_HINTS_METHOD} \
          ${GU_REFERENCE_DTSEQ} ${GU_REFERENCE_METHOD} \
          ${GU_INITIAL_OPT} ${GU_INITIAL_NOISE} \
          ${GU_MAXLEVEL}

clean::
	/bin/rm -fv out/*-*.{fni,txt,png}
