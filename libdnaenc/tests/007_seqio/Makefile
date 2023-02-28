# Last edited on 2023-02-27 19:21:41 by stolfi

TEST_LIB  := libdnaenc.a
TEST_LIB_DIR  := ../..

PROG := test_dnae_seqio

JS_LIBS := \
  libmsmatch.a \
  libgeo.a \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

# check:  clean run show
check:  clean run

SEQDIR := ${STOLFIHOME}/programs/c/DNA/DATA/small/bas

SEQ := 250A

INNAME := ${SEQDIR}/${SEQ}

OUTNAME := out/${SEQ}

clean::
	-/bin/rm -f ${OUTNAME}-*.pgm
  
run: ${PROG} out
	${PROG} \
          ${INNAME} \
          ${OUTNAME}

out:
	mkdir -p out

show: 
	@echo "** show what?"
