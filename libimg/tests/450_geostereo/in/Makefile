# Makefile to process POV-Ray lab hand-in files - STATIC IMAGES
# Last edited on 2020-10-11 15:30:02 by jstolfi

# Small square
# 
# WIDTH := 400
# HEIGHT := 400
# 
# ICONWIDTH := 80
# ICONHEIGHT := 80

# Large square

WIDTH := 640
HEIGHT := 640

ICONWIDTH := 80
ICONHEIGHT := 80

NRAYS = 3

OTHERINPUTS := \
  eixos.inc \
  camlight.inc

PAIR_NAME = A
POV_NAME = stereopair_${PAIR_NAME}

# No need to change below here

# POV-Ray installation directories:
POVRAY := /usr/bin/povray
POVINC := /usr/share/povray-3.7/include/

TOPDIR := ${STOLFIHOME}/public_html/courses/MC937-2017-1-A/progs
POVTTF := ${TOPDIR}/hand-out/tt-fonts

# These are defined with "=" rather than ":=" in order to defer expansion
DIRTAIL = ${shell pwd | sed -e 's:^.*/::'}
IMVIEW = display -title "${DIRTAIL}/%f"

.PHONY: all single cleanup  image showimage paramfile imagepair showpair 

all: check

check:  imagepair showpair

cleanup: 
	-/bin/rm -f ${PAIR_NAME}-*.{png,pgm,ppm}

imagepair:
	for eye in 0 1 ; do \
	  rm -f params.inc ; \
          ${MAKE} EYE=$${eye} single ; \
        done

showpair: ${PAIR_NAME}-0.png ${PAIR_NAME}-1.png
	${IMVIEW} `ls ${PAIR_NAME}-*.{png,pgm,ppm} | sort ` 

#ifneq "/${EYE}" "/"
#######################################################################
# When ${EYE} is defined

single: image showimage

IMG_NAME := ${PAIR_NAME}-${EYE}

image: ${IMG_NAME}.png ${IMG_NAME}.pgm ${IMG_NAME}.ppm 

showimage: ${IMG_NAME}.png
	-if [[ -s ${IMG_NAME}.png ]]; then ${IMVIEW} ${IMG_NAME}.png ; fi

${IMG_NAME}.png: ${POV_NAME}.pov ${OTHERINPUTS} paramfile
	-/bin/rm -f ${IMG_NAME}.png
	nice ${POVRAY} \
            +FN +Q9 \
            +W${WIDTH} +H${HEIGHT} \
            +AM1 +A0.0 +R${NRAYS} \
            +D \
            +L${POVINC} \
            +L${POVTTF} \
	    +I${POV_NAME}.pov \
	    +O${IMG_NAME}.png
        
paramfile: 
	echo "#declare eye = ${EYE};" > params.inc

${IMG_NAME}.ppm: ${IMG_NAME}.png
	convert  ${IMG_NAME}.png ${IMG_NAME}.ppm

${IMG_NAME}.pgm: ${IMG_NAME}.png
	convert  ${IMG_NAME}.png -colorspace Gray ${IMG_NAME}.pgm

#endif
