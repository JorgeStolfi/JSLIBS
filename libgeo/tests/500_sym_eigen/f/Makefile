# Last edited on 2024-11-22 20:53:33 by stolfi

PROG = test_sum_eigen_fortran

.PHONY:: all clean show

OUT_DIR := out
OUT_NAME := ${OUT_DIR}/${PROG}
  
all: check show

check:  ${OUT_NAME}.txt 

# DAMN MINUS ZEROS!!!!

${OUT_NAME}.txt: ${PROG} Makefile
	mkdir -p ${OUT_DIR}
	( ${PROG} 2>&1 ) | sed -e 's:[-]\(0[.]000*\b\): \1:g' > ${OUT_NAME}.txt

show: ${OUT_NAME}.txt
	@echo "=== ${OUT_NAME}.txt ==="
	cat ${OUT_NAME}.txt

clean::
	-/bin/rm -f ${PROG} *.o ${OUT_NAME}.txt

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
