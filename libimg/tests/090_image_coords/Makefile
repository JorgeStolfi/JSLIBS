# Last edited on 2024-12-05 07:58:20 by stolfi

PROG = test_image_coords

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run-tests 

all: check

check:  clean run-tests

run-tests: ${PROG}
	${PROG} \
          -xAxis right \
          -yAxis up \
          -iUnit 0.5 \
          -iOrg 10 30 \
          -oUnit 0.25 \
          -oCenter \
          -oSize 200 300 \
          -maxval 65535 \
          -verbose
          

        
clean::
