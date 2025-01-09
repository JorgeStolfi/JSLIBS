#! /bin/bash
# Last edited on 2025-01-07 22:20:56 by stolfi

nf="$1"; shift;

PROG="test_integrate_recursive"

printf "function ${nf} ..."
${PROG} -function ${nf}
for tag in IZ KZ  IG JG KG DG  IW JW KW DW ; do
  fni_view -scale 1.0 out/${nf}-*-*-${tag}.fni
done
echo "OK"
