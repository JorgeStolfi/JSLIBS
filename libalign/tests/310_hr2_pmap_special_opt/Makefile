# Last edited on 2024-11-08 11:18:15 by stolfi

PROG = test_hr2_pmap_special_opt

TEST_LIB := libalign.a
TEST_LIB_DIR := ../..

JS_LIBS :=  \
  libminn.a \
  libminu.a \
  libgeo.a \
  libimg.a \
  libjs.a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
.PHONY:: do-test
 
all: check

check:  
	for type in TRANSLATION CONGRUENCE SIMILARITY AFFINE ; do \
          ${MAKE} do-test TYPE=$$type ; \
        done

ifneq "/${TYPE}" "/"

OUT_PREFIX := out/test-${TYPE}

do-test: ${PROG}
	rm -f ${OUT_PREFIX}-*.{txt,dat,png}
	${PROG} ${TYPE}
	for datafile in ${OUT_PREFIX}-*-1D-plot.txt ; do \
          plot_opt_mismatch_1D.sh $$datafile ; \
        done
	display -title '%f' ${OUT_PREFIX}-*.png 

endif
