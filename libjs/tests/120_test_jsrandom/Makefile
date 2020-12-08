# Last edited on 2020-12-06 18:20:39 by jstolfi

PROG = test_jsrandom

TEST_LIB := libjs.a
TEST_LIB_DIR := ../..

JS_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test show-plots

all: check

check:  do-test show-plots

do-test: ${PROG}
	-rm -f out/*.{his,png}
	${PROG}

show-plots: out/drandom_qt_A.his
	plot_loggauss_hist.sh dloggaussrand_qt_A
	plot_unif_hist.sh drandom_qt_A
	plot_gauss_hist.sh dgaussrand_qt_A
	for xtype in int uint ; do \
          for size in 32 64 ; do \
            for end in lo hi ; do \
              plot_int_hist.sh $${xtype}$${size}_random_$${end}_A ; \
            done ; \
            for xtest in A B ; do \
              plot_int_hist.sh $${xtype}$${size}_abrandom_qt_$${xtest} ; \
            done \
          done \
        done

clean::
	rm -f out/*.his out/*.png
