#! /bin/bash 
# Last edited on 2025-02-28 14:21:52 by stolfi

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
  # in 1x1 2x2 4x3 16x12
  for size in 32x24 256x192 ; do
    nx="${size/x*/}"; ny="${size/*x/}"
    run_one_test.sh ${func_num} ${func_name} ${nx} ${ny} ${noisy} ${export}
  done
done
echo "OK"

