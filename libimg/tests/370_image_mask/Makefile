# Last edited on 2020-10-11 15:27:36 by jstolfi

PROG = test_image_mask

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test show-masks plot-masks plot-mask-single view-masks view-mask-single

all: check

# check:  clean do-test show-masks plot-masks
check:  clean do-test show-masks view-masks

OUTNAME := out/msk

do-test: ${PROG}
	${PROG} 400 400

show-masks:
	display ${OUTNAME}-*G.pgm
        
plot-masks:
	for f in out/msk-*L.pgm ; do \
	  ${MAKE} MASK=$$f plot-mask-single ; \
        done
        
view-masks:
	for f in out/msk-*L.pgm ; do \
	  ${MAKE} MASK=$$f view-mask-single ; \
        done
        
clean::
	/bin/rm -fv ${OUTNAME}-*.pgm

########################################################################        
# Secrtion for recursive "make" -- client must defined ${MASK}

ifneq "/${MASK}" "/"

TEMP := /tmp/tim

plot-mask-single: 
	pnm_to_fni \
            -min 0 -max 1 \
          < ${MASK} \
          > ${TEMP}.fni
	fni_plot.sh \
            -title "${MASK}" \
            -range -0.01 1.01 \
          < ${TEMP}.fni \
          > ${TEMP}.eps
	display -background white ${TEMP}.eps

view-mask-single: 
	fni_view \
            -isMask T \
            -colorize 0 1 \
            -scale 200 \
          ${MASK}

endif
# End ${MASK} section
########################################################################
