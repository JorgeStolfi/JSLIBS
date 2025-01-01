# Makefile for generic C program
# Last edited on 2024-12-31 20:54:38 by stolfi

PROG := test_cpk_plot

TEST_LIB = libcpk.a
TEST_LIB_DIR = ../..

JS_LIBS := \
  libeps.a \
  libgeo.a \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run

all: check

# check:  clean run
check:  run

run: ${PROG}
	run_tests.sh
