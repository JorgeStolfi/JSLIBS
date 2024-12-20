# Last edited on 2024-12-20 18:26:15 by stolfi

PROG = test_interpolate

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libgausol.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

NAME := ta
NC := 3
EXT := ppm

do-test: ${PROG} in/${NAME}-orig.${EXT}
	mkdir -p out
	${PROG} ${NAME} ${NC}
	display \
            -title '%f' \
            -filter Point -resize '800%' \
            out/${NAME}-orig.${EXT} \
            `ls out/${NAME}-is-*.${EXT} | sort` \
            `ls out/${NAME}-bd-*.${EXT} | sort`
	show-plot.sh out/${NAME}-RS-0003-0003.txt
	show-plot.sh out/${NAME}-RS-0000-0000.txt
	show-plot.sh out/${NAME}-RE-0000-0000.txt
	show-plot.sh out/${NAME}-RR-0000-0000.txt
	show-plot.sh out/${NAME}-RM-0000-0000.txt
	show-plot.sh out/${NAME}-RP-0000-0000.txt

clean::
	/bin/rm -fv out/*.ppm out/*.pgm out/*.txt out/*.fni
