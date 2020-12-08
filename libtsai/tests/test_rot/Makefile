# Last edited on 2020-10-11 15:33:27 by jstolfi

PROG := test_trans_matrix
TEST_LIB := libtsai.a
TEST_LIBDIR := ../..

NON_PROGS :=

IGNORE :=  

JS_LIBS:= \
  liblmdif.a \
  libgeo.a \
  libjs.a

all: check

check:  build run_b

# OTHER_C_FLAGS := -pg
# OTHER_LD_FLAGS := -pg

include ../../../GENERIC-LIB-TEST.make

run_a: ${PROG} 
	rm -rf out/* 
	./${PROG} \
	../../../data/videos/stage/povray/stage02/out/out \
	svga \
	../../../programs/camera_gen/stage02/gcalibration.dat \
        ../../../programs/affine_tracker/data/stage02/00-18-world-coords.dat	

	#../../../programs/camera_gen/stage01/gcalibration.dat 
	#fcalibration_stage.txt 
run_b: ${PROG} 
	rm -rf out/* 
	./${PROG} \
	../../../data/videos/stage/povray/stage01/out/out \
	hvga \
	../../../programs/affine_tracker/fcalibration_sem_kappa.txt \
	../../../programs/affine_tracker/data/stage01/00-23-world-coords.dat 	

run_c: ${PROG} 
	rm -rf out/* 
	./${PROG} \
	../../../data/videos/office/modified_office/out/out \
	hvga \
	../../../programs/camera_gen/office01/gcalibration.dat \
	../../../programs/affine_tracker/data/office01/00-23-world-coords.dat 	

bug: ${PROG} 
	rm -rf out/* 
	./${PROG} \
	../../../data/videos/stage/real/0149/out \
	optura \
	 ../../../programs/affine_tracker/fcalibration.txt \
	../../../programs/affine_tracker/data/0149/00-23-world-coords.dat 	

bug2: ${PROG}
	rm -rf out/* 
	./${PROG} \
	../../../data/videos/stage/real/0161/out \
	optura \
	../../../programs/affine_tracker/fcalibration.txt \
	../../../programs/affine_tracker/data/0161/00-23-world-coords.dat 	


bug1: ${PROG} 
	./${PROG} \
	../../../data/videos/stage/povray/stage01/out/out \
	optura \
	../../../programs/camera_gen/stage01/gcalibration.dat  \
	../../../programs/affine_tracker/data/0152/00-23-world-coords.dat 	

#_artoolkit.txt
#../../../programs/PaperComparision/gcalibration.dat 

stage04: ${PROG} 
	./${PROG} \
	../../../data/videos/stage/povray/stage04/out_all/out \
	svga \
	../../../programs/affine_tracker/fcalibration.txt \
	../../../programs/affine_tracker/data/stage04/00-08-world-coords-artoolkit.dat 	

artoolkit: ${PROG} 
	./${PROG} \
	../../../data/videos/stage/povray/stage05/out/out \
	svga \
	../../../programs/affine_tracker/fcalibration.txt \
	../../../programs/affine_tracker/data/stage04/00-08-world-coords-artoolkit.dat 	


#Na ultima linha passa o cpar real ou none quando este nao existe
office: ${PROG}
	rm -rf points
	mkdir points
	./${PROG} \
	../../../data/videos/office/modified_office/OFFICE/out \
	svga \
	/home2/minetto/affine_tracker/programs/Softwares/Tsai/Tsai-method-v3.0b3/calibration_tsai.txt \
	../../../programs/affine_tracker/data/office/world-coords.dat \
	../../../programs/camera_gen/office01/gcalibration.dat 

n8: ${PROG}
	rm -rf out/*
	rm -rf points
	mkdir points
	./${PROG} \
	../../../data/videos/stage/real/0149/out \
	optura \
	in/t8/opt.cpar \
	in/t8/p_w.txt \
	none


