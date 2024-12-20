# Last edited on 2024-12-08 18:51:06 by stolfi

SHELL := /bin/bash

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_raytrace

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
