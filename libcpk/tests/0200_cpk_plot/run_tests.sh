#! /bin/bash
# Last edited on 2024-12-31 20:58:56 by stolfi

# Quais conjuntos de dados devem ser usados nos testes:
datasets=( square seila2 seila1 naopau saopau almira )

for dataset in ${datasets[@]}; do
  mkdir -p out
  ( cd out && rm -f ${dataset}.{eps,sol} )
  cmd=( test_cpk_plot ${options[@]} data/${dataset} out/${dataset} )
  echo "${cmd[@]}"
  ${cmd[@]}
  for tag in nd wd; do 
    if [ -r out/${dataset}/greed-${tag}.eps ]; then
      evince out/${dataset}/greed-${tag}.eps &
      evince out/${dataset}/grasp-${tag}.eps
    else
      echo "file out/${dataset}/greed-${tag}.eps not found"
    fi
  done
done
