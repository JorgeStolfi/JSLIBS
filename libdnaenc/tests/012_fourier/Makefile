# Last edited on 2023-02-14 07:30:10 by stolfi

TEST_LIB := libdnaenc.a
TEST_LIB_DIR := ../..
PROG := test_dnae_fourier

JS_LIBS := \
  libmsmatch.a \
  libgeo.a \
  libeps.a \
  libimg.a \
  libjs.a

OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libfftw3.so

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: run show-main show-all depend

all: check

check:  clean run show-main
# check:  clean run

DATADIR := ../00-DATA

MAX_LEVEL := 5
NUM_POINTS := 2048

INIT_WEIGHTS := wg020
INCR_WEIGHTS := wg060

PLOT_WIDTH := 72.0
PLOT_HEIGHT := 18.0
FONT_SIZE := 10.0

OUTNAME := out/out

SHOWPS := evince

clean:: 
	-/bin/rm -f ${OUTNAME}-*.{pgm,egs,eps}
  
run: ${PROG} 
	${PROG} \
          -maxLevel ${MAX_LEVEL} \
          -numPoints ${NUM_POINTS} \
          -initFilter `cat ${DATADIR}/weights/${INIT_WEIGHTS}.txt | sed -e 's:[\#].*$$::g'` \
          -incrFilter `cat ${DATADIR}/weights/${INCR_WEIGHTS}.txt | sed -e 's:[\#].*$$::g'` \
          -plotSize ${PLOT_WIDTH} ${PLOT_HEIGHT} \
          -fontSize ${FONT_SIZE} \
          ${OUTNAME}

show-main: ${OUTNAME}-00-db.eps
	@for psf in \
          ` ls ${OUTNAME}-??-fo.eps ${OUTNAME}-??-fr.eps ${OUTNAME}-??-fs.eps | sort ` ; \
        do \
          echo "$$psf" ; ${SHOWPS} $$psf ; \
        done

show-all:  ${OUTNAME}-00-db.eps
	@for psf in \
          ` ( ls ${OUTNAME}-??-{db,f0,f1}.eps | sort ) ; \
            ( ls ${OUTNAME}-??-fo.eps ${OUTNAME}-??-fr.eps ${OUTNAME}-??-fs.eps | sort ) \
          ` ; do \
          echo "$$psf" ; ${SHOWPS} $$psf ; \
        done
