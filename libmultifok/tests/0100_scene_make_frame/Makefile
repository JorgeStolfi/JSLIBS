# Makefile para exercícios ESTÁTICOS de POV-Ray
# Last edited on 2024-12-16 18:10:22 by stolfi

SHELL := /bin/bash

TEST_LIB := libmultifok.a
TEST_LIB_DIR := ../..
PROG := test_mfok_scene_make_frame

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

all: make_images.py
	make_images.py

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
