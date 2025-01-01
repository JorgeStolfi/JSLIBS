# Makefile for generic C program
# Last edited on 2025-01-01 02:06:01 by stolfi

PROG := test_hxg

TEST_LIB = libhxg.a
TEST_LIB_DIR = ../..

JS_LIBS := \
  libeps.a \
  libgeo.a \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run

all: run 

run: ${PROG} 
	${PROG}
	evince out/test.eps

