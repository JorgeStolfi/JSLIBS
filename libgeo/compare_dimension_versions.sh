#! /bin/bash 
# Last edited on 2024-08-30 11:43:59 by stolfi

# Compares each of "r2.h", "r3.h", "r2x2.h", "r3x3.h" "hr2.h" "hr3.h"
# with the version of the same file, one dimension up,
# after doing the obvious substitutions on the former.
# E.g. compares "r2.h" with "r3.h", "r3x3.h" with "r4x4.h", etc.
#
# The comparison results are written to ".{ff}.diff" where
# {gg} is the second module.

tmp="/tmp/$$"

for fg in \
    r2.h:r3.h r2x2.h:r3x3.h hr2.h:hr3.h \
    r3.h:r4.h r3x3.h:r4x4.h hr2.h:hr3.h \
  ; do 
  ff="${fg/:*/}" 
  gg="${fg/*:/}" 
  cat ${ff} \
    | sed -f up_one_dimension.sed \
    > ${tmp}-${ff}.upped
  prdiff -Bb ${tmp}-${ff}.upped ${gg} \
    > .${gg}.diff
done 

rm ${tmp}-*
