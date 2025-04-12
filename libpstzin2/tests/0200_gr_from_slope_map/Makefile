# Last edited on 2025-03-15 04:43:52 by stolfi

PROG = test_gr_from_slope_map

TEST_LIB := libpstzin2.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libpst.a \
  libeps.a \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-tests

all: ${PROG}
	${PROG} -function 17 cbabel -noisy N -size 64 48 
	evince out/17*-graph.eps

# all: check
# 
# check:  clean do-tests
# 
# ZFUNCS := 01 02 03 04 05 11 
# 
# do-tests: ${PROG}
# 	mkdir -p out
# 	for nf in ${ZFUNCS} ; do \
# 	  ./run_test.sh $$nf ; \
#         done

clean::
	/bin/rm -fv out/*-*.fni
