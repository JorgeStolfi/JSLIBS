# Last edited on 2023-12-18 10:57:33 by stolfi

TEST_LIB := libneuro.a
TEST_LIB_DIR := ../..
PROG := test_filter

JS_LIBS := \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libfftw3.so \

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

# Nominal sampling freq (Hz): 
FSMP := 1.000 

# Nominal cutoff and shoulder freqs for bandpass:
FLO0 := 0.0030
FLO1 := 0.0120
FHI1 := 0.1200
FHI0 := 0.4800

check:  ${PROG}
	mkdir -p out
	rm -f out/*.txt out/*.png
	./${PROG} ${FSMP} ${FLO0} ${FLO1} ${FHI1} ${FHI0}
	for f in out/*_tab_*.txt ; do \
          ./plot_tabulated_filter.sh $$f ${FSMP} ; \
        done
	for f in out/*{lopa,band}.txt ; do \
	  ./plot_filter.sh $$f ${FSMP} ${FLO0} ${FLO1} ${FHI1} ${FHI0} ; \
        done
