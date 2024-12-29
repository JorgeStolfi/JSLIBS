#! /bin/bash
# Last edited on 2024-12-23 09:16:22 by stolfi

nf="$1"; shift;

PROG="test_maps_expand_reduce"

printf "function ${nf} ..."
${PROG} -function ${nf}
for tag in IZ KZ  IG JG KG DG  IW JW KW DW ; do
  fni_view -scale 1.0 out/${nf}-*-*-${tag}.fni
done
echo "OK"
