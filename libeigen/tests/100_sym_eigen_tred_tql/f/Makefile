# Last edited on 2024-12-05 17:12:52 by stolfi

SHELL := /bin/bash
PROG = test_sym_eigen_tred_tql_F

.PHONY:: all clean show

OUT_DIR := out
OUT_FILE_1 := ${OUT_DIR}/version1.txt
OUT_FILE_2 := ${OUT_DIR}/version2.txt
  
all: check

${OUT_FILE_1}: check
${OUT_FILE_2}: check

check: ${PROG} Makefile
	mkdir -p ${OUT_DIR}
	rm -f ${OUT_FILE_1} ${OUT_FILE_2}
	${PROG} > out/.log 2>&1 
	@for ff in ${OUT_FILE_1} ${OUT_FILE_2} ; do \
	  gg="$${ff/out/ref}" ; \
	  echo "=== $${ff} $${gg} ===" ; \
          prdiff -Bb $${ff} $${gg} ; \
        done

show: ${OUT_FILE_1} ${OUT_FILE_2}
	@for ff in ${OUT_FILE_1} ${OUT_FILE_2} ; do \
	  echo "=== $${ff} ===" ; \
	  cat $${ff} ; \
        done

clean::
	-/bin/rm -f ${PROG} *.o ${OUT_FILE_1} ${OUT_FILE_2}

# ======================================================================
# Rules to compile the FORTRAN sources

FFILES := tred1.f tred2.f tred1m.f tred2m.f tql1.f tql1m.f tql2.f tql2m.f pythag.f

OFILES := ${subst .f,.o,${FFILES}}

LIBS :=
  
# FORTRAN 77 compiler:
FC := gfortran -std=legacy -fdebug-aux-vars -g -O0

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
