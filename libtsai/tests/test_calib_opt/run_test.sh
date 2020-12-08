#! /bin/bash
# Last edited on 2015-10-04 20:43:57 by stolfilocal

tdir="$1"; shift; # Test data subdirectory "t01", "a01", etc.

PROG=test_calib_opt
PSVIEW=okular

mkdir -p out/${tdir}
rm -rf out/${tdir}/*
./${PROG} \
    `cat in/${tdir}/model.txt` \
    in/${tdir} \
    out/${tdir} 
    
${PSVIEW} out/${tdir}/true_guess.eps
${PSVIEW} out/${tdir}/true_opt1.eps
${PSVIEW} out/${tdir}/true_opt2.eps
