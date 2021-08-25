# Last edited on 2021-08-25 03:23:24 by stolfi

PROG := test_neuromat_eeg_image

TEST_LIB := libneuro.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib/x86_64-linux-gnu/libpng.so \
  /usr/lib/x86_64-linux-gnu/libz.so

all: check

check:  ${PROG}
	rm -f out/*.{fni,png}
	${PROG}
	cd out && for f in *.fni ; do \
          fni_to_png < $$f > $${f/.fni/.png} ; \
        done ; \
        eom *.png

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

clean:: 
	rm -f out/*.{fni,png}
