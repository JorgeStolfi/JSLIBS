# Last edited on 2024-01-06 08:28:04 by stolfi

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
FLO0 := 0.0050
FLO1 := 0.0150
FHI1 := 0.1600
FHI0 := 0.4800

check:  ${PROG}
	mkdir -p out
	rm -f out/*.txt out/*.png
	./${PROG} ${FSMP} ${FLO0} ${FLO1} ${FHI1} ${FHI0}
	@for ftype in LG ER GA BW BQ ; do \
          for band in lopa band ; do \
            for f in out/$${ftype}_$${band}.txt ; do \
	      if [[ -s $${f} ]]; then \
                ./plot_filter.sh 0 $$f ${FSMP} ${FLO0} ${FLO1} ${FHI1} ${FHI0} ; \
                ./plot_filter.sh 1 $$f ${FSMP} ${FLO0} ${FLO1} ${FHI1} ${FHI0} ; \
              fi ; \
            done ; \
          done ; \
        done
	for f in out/*_tab_*.txt ; do \
          ./plot_tabulated_filter.sh $$f ${FSMP} ; \
        done
	./plot_clear_tiny.sh out/clear_tiny.txt 2.0
	./plot_clear_tiny.sh out/clear_tiny.txt 1.0e-8
