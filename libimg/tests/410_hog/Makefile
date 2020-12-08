# Last edited on 2020-10-11 15:40:54 by jstolfi

PROG = test_hog

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

TESTS := \
  bdots/036/0.010 \
  bpent/045/0.010 \
  btria/036/0.010 \
  bueye/036/0.010 \
  gras2/036/0.010 \
  grass/036/0.010 \
  janus/036/0.010 \
  ptpyr/045/0.010 \
  wavec/036/0.010 \

.PHONY:: all run-all run-single clean-all clean-single plot-approx

all: check

check:  run-all
clean:: clean-all

run-all:  
	for t in ${TESTS} ; do \
          ${MAKE} TEST="$${t}" run-single ; \
        done

clean-all: 
	for t in ${TESTS} ; do \
          ${MAKE} TEST="$${t}" clean-single ; \
        done


TEST := 
ifneq "/${TEST}" "/"
#######################################################################
# Recursive "make" to run a single test

TESTARGS := ${subst /, ,${TEST}}

PREFIX         := ${word 1, ${TESTARGS}}
NBINS          := ${word 2, ${TESTARGS}}
NOISE          := ${word 3, ${TESTARGS}}
 
RDFILES := \
  in/${PREFIX}-img.ppm \
  in/${PREFIX}-msk.pgm
  
GFILE_FNI := out/${PREFIX}-G.fni
GFILE_PPM := out/${PREFIX}-G.ppm
  
SHOW := display -title '%f'
        
run-single: ${PROG} ${RDFILES}
	-rm -f out/${PREFIX}-*.*
	${PROG} ${PREFIX} ${NOISE} ${NBINS}
	plot-gradient.sh ${GFILE_FNI}
	fni_to_pnm -min 0 -1 -1 -max 1 +1 +1 < ${GFILE_FNI} > ${GFILE_PPM}
	fni_view ${GFILE_FNI}
	for ff in out/${PREFIX}-*.txt ; do \
          echo "=== $${ff} ====================" ; \
          cat $${ff} ; \
	  plot-hog.sh ${NBINS} $${ff} ; \
        done

clean-single::
	/bin/rm -f out/${PREFIX}-*.txt out/${PREFIX}-*.png out/${PREFIX}-*.fni out/${PREFIX}-*.ppm
 
# 
# End make with ${TEST}
######################################################################
endif
       
# ======================================================================        
# Trying to find a formula for the azimuth distribution.

NA := 360

plot-approx: compute-polar-gaussian.gawk plot-polar-gaussian.sh
	compute-polar-gaussian.gawk -v na="${NA}" > out/plot.data
	plot-polar-gaussian.sh ${NA} out/plot.data
