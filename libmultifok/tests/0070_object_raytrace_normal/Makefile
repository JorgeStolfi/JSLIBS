# Makefile para exercícios ESTÁTICOS de POV-Ray
# Last edited on 2025-02-03 02:31:04 by stolfi

SHELL := /bin/bash

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..

PROG := test_mfok_object_raytrace_normal

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libgausol.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libjpeg.so \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so
  
# Profiling:
OTHER_C_FLAGS := -pg
OTHER_LD_FLAGS := -pg

all: ${PROG}
	${PROG}

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
