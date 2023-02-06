# Last edited on 2023-02-04 07:05:40 by stolfi

PROG := test_calib_guess
TEST_LIB := libtsai.a
TEST_LIB_DIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS:= \
  libgeo.a \
  libeps.a \
  libjs.a
  
all: check

check:  build run-tests

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

# PSVIEW := evince
PSVIEW := evince

include ../../../../GENERIC-LIB-TEST.make

run-tests: t03-1
# BUG: run-tests: t03-2 

# Tests ----------------------------------------------------------------------

t03-1: ${PROG}
	-rm -rf out/t03-1/*
	./${PROG} \
	    svga \
	    in/t03 \
            out/t03-1 \
            1 
	${PSVIEW} out/t03-1/true_affine.eps out/t03-1/true_guess.eps

t03-2: ${PROG}
	-rm -rf out/t03-2/*
	./${PROG} \
	    svga \
	    in/t03 \
            out/t03-2 \
            2 
	${PSVIEW} out/t03-2/true_affine.eps out/t03-2/true_guess.eps

n1: ${PROG}
	-rm -rf out/n1/*
	./${PROG} \
	    optura \
	    in/t1 \
            out/n1 \
            2  
	${PSVIEW} out/n1/true_affine.eps
	${PSVIEW} out/n1/true_guess.eps

c2: ${PROG}
	-rm -rf out/c2/*
	./${PROG} \
	    optura \
	    in/t2 \
            out/c2 \
            2  
	${PSVIEW} out/c2/true_affine.eps
	${PSVIEW} out/c2/true_guess.eps

c4: ${PROG}
	-rm -rf out/*
	./${PROG} \
	    svga \
	    in/t4 \
            out/c4 \
            2  
	${PSVIEW} out/c4/true_affine.eps
	${PSVIEW} out/c4/true_guess.eps

n5: ${PROG}
	-rm -rf out/n5/*
	./${PROG} \
	    svga \
	    in/t5 \
            out/n5 \
            2  
	${PSVIEW} out/n5/true_affine.eps
	${PSVIEW} out/n5/true_guess.eps

n6: ${PROG}
	-rm -rf out/n6/*
	./${PROG} \
	    svga \
	    in/t6 \
            out/n6 \
            2  
	${PSVIEW} out/n6/true_affine.eps
	${PSVIEW} out/n6/true_guess.eps

n7: ${PROG}
	-rm -rf out/n7/*
	./${PROG} \
	    optura \
	    in/t7 \
            out/n7 \
            2  
	${PSVIEW} out/n7/true_affine.eps
	${PSVIEW} out/n7/true_guess.eps

n8: ${PROG}
	-rm -rf out/n8/*
	./${PROG} \
	    optura \
	    in/t8 \
            out/n8 \
            2  
	${PSVIEW} out/n8/true_affine.eps
	${PSVIEW} out/n8/true_guess.eps

n9: ${PROG}
	-rm -rf out/n9/*
	./${PROG} \
	    hvga \
	    in/t9 \
            out/n9 \
            2  
	${PSVIEW} out/n9/true_affine.eps
	${PSVIEW} out/n9/true_guess.eps

n10: ${PROG}
	-rm -rf out/n10/*
	./${PROG} \
	    hvga \
	    in/t10 \
            out/n10 \
            2  
	${PSVIEW} out/n10/true_affine.eps
	${PSVIEW} out/n10/true_guess.eps

n11: ${PROG}
	-rm -rf out/n11/*
	./${PROG} \
	    hvga \
	    in/t11 \
            out/n11 \
            2  
	${PSVIEW} out/n11/true_affine.eps
	${PSVIEW} out/n11/true_guess.eps

n12: ${PROG}
	-rm -rf out/n12/*
	./${PROG} \
	    hvga \
	    in/t12 \
            out/n12 \
            2  
	${PSVIEW} out/n12/true_affine.eps
	${PSVIEW} out/n12/true_guess.eps

n13: ${PROG}
	-rm -rf out/n13/*
	./${PROG} \
	    sony \
	    in/t13 \
            out/n13 \
            2  
	${PSVIEW} out/n13/true_affine.eps
	${PSVIEW} out/n13/true_guess.eps

n14: ${PROG}
	-rm -rf out/n14/*
	./${PROG} \
	    sony \
	    in/t14 \
            out/n14 \
            2  
	${PSVIEW} out/n14/true_affine.eps
	${PSVIEW} out/n14/true_guess.eps


a2: ${PROG}
	-rm -rf out/a2/*
	./${PROG} \
	    svga \
	    in/a2 \
            out/a2 \
            2  
	${PSVIEW} out/a2/true_affine.eps
	${PSVIEW} out/a2/true_guess.eps

a3: ${PROG}
	-rm -rf out/a3/*
	./${PROG} \
	    hvga \
	    in/a3 \
            out/a3 \
            2  
	${PSVIEW} out/a3/true_affine.eps
	${PSVIEW} out/a3/true_guess.eps

a4: ${PROG}
	-rm -rf out/a4/*
	./${PROG} \
	    hvga \
	    in/a4 \
            out/a4 \
            2  
	${PSVIEW} out/a4/true_affine.eps
	${PSVIEW} out/a4/true_guess.eps


a5: ${PROG}
	-rm -rf out/a5/*
	./${PROG} \
	    hvga \
	    in/a5 \
            out/a5 \
            2  
	${PSVIEW} out/a5/true_affine.eps
	${PSVIEW} out/a5/true_guess.eps

a6: ${PROG}
	-rm -rf out/a6/*
	./${PROG} \
	    hvga \
	    in/a6 \
            out/a6 \
            2  
	${PSVIEW} out/a6/true_affine.eps
	${PSVIEW} out/a6/true_guess.eps


