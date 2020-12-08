# Last edited on 2020-10-11 15:38:10 by jstolfi

PROG = testf

.PHONY:: all clean

OUT_DIR := out
OUT_NAME := ${OUT_DIR}/${PROG}
  
all: check

check:  ${OUT_NAME}.txt 

${OUT_NAME}.txt: ${PROG} Makefile
	mkdir -p ${OUT_DIR}
	${PROG} > ${OUT_NAME}.txt 2>&1
	@echo "=== ${OUT_NAME}.txt ==="
	cat ${OUT_NAME}.txt

clean::
	-/bin/rm -f ${OUT_NAME}.txt

# ======================================================================
# Rules to compile the FORTRAN sources

FFILES := tred1.f tred2.f tql1.f tql2.f pythag.f

OFILES := ${subst .f,.o,${FFILES}}

LIBS :=
  
# FORTRAN 77 compiler:
FC := gfortran -std=legacy

%.o: %.f
	${FC} -c -o $*.o $*.f 

build: ${PROG}

${PROG}: ${PROG}.o ${OFILES} ${LIBS} Makefile
	${FC} -o ${PROG} ${PROG}.o ${OFILES} ${LIBS}


# ====================================================================
# Dependencies (by hand only):

depend: 
	@echo "no automatic dependencies"

testf.o:   testf.f
pythag.o:  pythag.f
tql1.o:    tql1.f
tql2.o:    tql2.f
tred1.o:   tred1.f
tred2.o:   tred2.f 
