# Last edited on 2025-01-16 06:13:38 by stolfi

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

check:  clean do-tests

# TESTS := 01-ramp10-N 03-ramp11-N 
# TESTS := 03-ramp11-N 
TESTS := 04-parabo-N

SIZE := 0016x0012
# SIZE := 0256x0192

do-tests: ${PROG}
	mkdir -p out
	for test in ${TESTS} ; do \
	  func="$${test}"; \
	  ./run_test.sh $${func} ${SIZE}; \
        done
 
clean::
	/bin/rm -fv out/*-*.{fni,txt,png}
