# Last edited on 2025-02-04 01:48:19 by stolfi

SHELL := /bin/bash

PROG := test_mfok_render

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

all: run_test.py
	run_test.py

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
