#! /bin/bash 
# Last edited on 2025-01-26 10:06:54 by stolfi

for func in `cat in/funcs.txt | sed -e 's:[#].*$::g'`; do 
  func_num="${func/-*/}"
  func_name="${func/*-/}"; func_name="${func_name/:*/}"
  noisy="${func/*:/}"
  #  1x1 2x2 4x3 16x12 32x24 
  for size in 256x192 ; do
    nx="${size/x*/}"; ny="${size/*x/}"
    run_one_test.sh ${func_num} ${func_name} ${nx} ${ny} ${noisy}
  done
done
echo "OK"

