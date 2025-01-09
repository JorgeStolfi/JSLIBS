# Last edited on 2025-01-09 05:49:42 by stolfi

PROG = test_integrate_iterative

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

check:  clean do-tests

# TESTS := 01-ramp10:F 03-ramp11:F 
# TESTS := 03-ramp11:F 
TESTS := 04-parabo:F

# SIZE := 0001x0001
SIZE := 0016x0012

do-tests: ${PROG}
	mkdir -p out
	for test in ${TESTS} ; do \
	  func="$${test/:*/}"; keepNull="$${test/*:/}" ; \
	  ./run_test.sh $${func} $${keepNull} ${SIZE}; \
        done
 
clean::
	/bin/rm -fv out/*-*.fni
