# Last edited on 2014-05-23 13:43:41 by stolfilocal

ALL_TESTS := \
  ellipse \
  r2 \
  r3 \
  r4 \
  r6 \
  rn \
  hr2 \
  hr3 \
  sym_eigen \
  gauss_elim \
  i2 \
  i3 \
  hi2
  
SOME_TESTS := \
  r2 \
  r3 \
  r4 \
  r6 \
  rn

SUBDIRS := ${ALL_TESTS}  

all: check

check:  build-tests check
  
include ${STOLFIHOME}/programs/c/GENERIC-ROOT-DIR.make
