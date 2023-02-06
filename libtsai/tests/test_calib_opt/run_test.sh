#! /bin/bash
# Last edited on 2023-02-04 08:15:54 by stolfi

tdir="$1"; shift; # Test data subdirectory "t01", "a01", etc.

PROG=test_calib_opt
PSVIEW=evince

mkdir -p out/${tdir}
rm -rf out/${tdir}/*
./${PROG} \
    `cat in/${tdir}/model.txt` \
    in/${tdir} \
    out/${tdir} 
    
outFile="out/${tdir}/true_guess.eps"
if [[ -s ${outFile} ]]; then
  ${PSVIEW} ${outFile}
  ${PSVIEW} out/${tdir}/true_opt1.eps
  ${PSVIEW} out/${tdir}/true_opt2.eps
else
  echo "** ${PROG} failed - ${outFile} not generated" 1>&2 ; exit 1
fi
