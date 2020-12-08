#! /bin/bash
# Last edited on 2011-07-07 22:10:32 by stolfilocal

nf="$1"; shift;

PROG="test_expand_reduce_maps"

printf "function ${nf} ..."
${PROG} -function ${nf}
for tag in IZ KZ  IG JG KG DG  IW JW KW DW ; do
  fni_view -scale 1.0 out/${nf}-*-*-${tag}.fni
done
echo "OK"
