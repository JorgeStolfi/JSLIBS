#! /bin/bash
# Last edited on 2024-01-08 17:04:41 by stolfi

for dfile_ini in out/train_*_ini_*.txt ; do
  echo "============================================" 1>&2
  dfile_fin="${dfile_ini/ini/fin}" 
  ls -l  ${dfile_ini} ${dfile_fin} 
  np="${dfile_ini/*train_/}" 
  np="${np/_*.txt/}" 
  for dfile in ${dfile_ini} ${dfile_fin} ; do
    plot_polygauss_train.sh NO ${dfile} ${np}
  done
  for ysc in 0 1 ; do
    pfile_ini="${dfile_ini/.txt/}_${ysc}.png" 
    pfile_fin="${dfile_fin/.txt/}_${ysc}.png" 
    ls -l  ${pfile_ini} ${pfile_fin} 
    if [[ ( -s ${pfile_ini} ) && ( -s ${pfile_fin} ) ]]; then
      pfile="${pfile_ini/_ini/}" 
      convert +append ${pfile_ini} ${pfile_fin} ${pfile} 
      display -title "%f" ${pfile} 
    fi
  done
done
