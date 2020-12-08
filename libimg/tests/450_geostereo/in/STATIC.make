# Makefile to process POV-Ray lab hand-in files - STATIC IMAGES
# Last edited on 2017-04-25 02:25:43 by stolfilocal

# Large portrait
# 
# WIDTH := 480
# HEIGHT := 640
# 
# ICONWIDTH := 75
# ICONHEIGHT := 100

# Small landscape
# 
# WIDTH := 320
# HEIGHT := 240
# 
# ICONWIDTH := 80
# ICONHEIGHT := 60

# Medium landscape
# 
# WIDTH := 480
# HEIGHT := 360
# 
# ICONWIDTH := 80
# ICONHEIGHT := 60

# Large landscape

WIDTH := 640
HEIGHT := 480

ICONWIDTH := 80
ICONHEIGHT := 60

NRAYS = 2

OTHERINPUTS := \
  camlight.inc

# No need to change below here

NAME = main

# POV-Ray installation directories:
POVRAY := /usr/bin/povray
POVINC := /usr/share/povray-3.7/include/

TOPDIR := ${STOLFIHOME}/public_html/courses/MC937-2017-1-A/progs
POVTTF := ${TOPDIR}/hand-out/tt-fonts

# These are defined with "=" rather than ":=" in order to defer expansion
DIRTAIL = ${shell pwd | sed -e 's:^.*/::'}
IMVIEW = display -title "${DIRTAIL}/%f"

.PHONY: all cleanup  image showimage all-images index index-images

all: image
# all: image showimage

cleanup: 
	-/bin/rm -f ${NAME}.png 
	-/bin/rm -f ${NAME}-*-i.png ${NAME}-i.png
	-/bin/rm -f ${NAME}-*.log ${NAME}.log

image: ${NAME}.png 

showimage: ${NAME}.png
	-if [[ -s ${NAME}.png ]]; then ${IMVIEW} ${NAME}.png ; fi

${NAME}.png: ${NAME}.pov ${OTHERINPUTS}
	-/bin/rm -f ${NAME}.png
	nice ${POVRAY} \
            +FN +Q9 \
            +W${WIDTH} +H${HEIGHT} \
            +AM1 +A0.0 +R${NRAYS} \
            +D \
            +L${POVINC} \
            +L${POVTTF} \
	    +I${NAME}.pov \
	    +O${NAME}.png
        
all-images: rausrs.txt
	../make-all-images.sh \
          `../list-ras.sh`

INDEXOPTS := -iconHeight ${ICONHEIGHT} 

INDEX_IMAGES := \
  absent.png absent-i.png \
  missing.png missing-i.png

index: title.txt rausrs.txt aus.txt ${INDEX_IMAGES}
	dir=`pwd`; \
        ../make-labday-image-page.sh ${INDEXOPTS} "`cat title.txt`"

index-orig: title.txt rausrs.txt aus.txt ${INDEX_IMAGES}
	../make-labday-image-page.sh ${INDEXOPTS} -absentOK -original "`cat title.txt`"

absent.png: Makefile
	ppmmake 'rgbi:0.75/0.75/0.75' ${WIDTH} ${HEIGHT} \
          | pnmdepth 255 \
          | ppmquant 254 \
          | convert PPM:- absent.png
	-${IMVIEW} absent.png
        
missing.png: Makefile
	ppmmake 'rgbi:0.80/0.15/0.15' ${WIDTH} ${HEIGHT} \
          | pnmdepth 255 \
          | ppmquant 254 \
          | convert PPM:- missing.png
	-${IMVIEW} missing.png

%-i.png: %.png
	convert $*.png \
            -geometry "x${ICONHEIGHT}" \
            +profile '*' \
            -colors 254 \
          $*-i.png
