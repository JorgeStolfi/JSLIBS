# Last edited on 2023-02-26 05:40:38 by stolfi

PROG = test_geostereo

TEST_LIB := libimg.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

.PHONY:: do-test

all: check

check:  clean do-test

NAME0 := A-0
NAME1 := A-1
NC := 3
EXT := ppm

do-test: ${PROG} in/${NAME0}.${EXT} in/${NAME1}.${EXT}
	mkdir -p out
	rm -f out/{disp,score}-*.*
	${PROG} ${NAME0} ${NAME1}
	display \
            -title '%f' \
            -filter Point -resize '200%' \
            in/${NAME0}.${EXT} \
            in/${NAME1}.${EXT} \
            out/disp-*.p?m \
            out/score-*.p?m

clean::
	/bin/rm -fv out/*.ppm out/*.pgm out/*.txt out/*.fni
