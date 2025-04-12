#! /bin/bash 
# Last edited on 2025-04-08 09:23:07 by stolfi

export="$1"; shift

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 1>&2
PROG="make_test_maps"
make -k -f Makefile ${PROG}
( cd ~/programs/c/C-IMG/fni_view   && make -k -f Makefile all )
( cd ~/programs/c/C-IMG/fni_to_pnm && make -k -f Makefile all )
( cd ~/programs/c/C-IMG/fni_hist   && make -k -f Makefile all )
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 1>&2

for func in `cat in/funcs.txt | sed -e 's:[#].*$::g'`; do 
  func_num="${func/-*/}"
  func_name="${func/*-/}"; func_name="${func_name/:*/}"
  noisy="${func/*:/}"
  # in 1x1 2x2 4x3 16x12 32x24 64x48 128x96 256x192 512x384
  for size in 16x12  ; do
    nx="${size/x*/}"; ny="${size/*x/}"
    run_one_test.sh ${func_num} ${func_name} ${nx} ${ny} ${noisy} ${export}
  done
done
echo "OK"

