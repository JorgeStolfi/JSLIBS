#! /bin/bash
# Last edited on 2024-12-31 17:41:13 by stolfi

# Quais conjuntos de dados devem ser usados nos testes:
datasets=( sq nosp sp a c altamira )
# datasets=( c )
# datasets=( altamira )
# datasets=( sp )
# datasets=( c sp )
# datasets=( sq )

# options=( -plot )
# options=( -plot -validate -mag4 )
# options=( -plot -validate -verbose -mag4 )
options=( -plot -validate )
# options=( -plot -validate -verbose )
# options=( -plot -validate -mag4 )
# options=( )

for dataset in ${datasets[@]}; do
  /bin/rm -f out/${dataset}/*.{ps,sol}
  cmd=( test_cpk ${options[@]} data/${dataset} out/${dataset} )
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
