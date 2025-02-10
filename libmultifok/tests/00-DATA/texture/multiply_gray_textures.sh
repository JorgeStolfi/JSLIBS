#! /bin/bash
# Last edited on 2025-02-08 20:18:26 by stolfi

file1="$1"; shift;
file2="$1"; shift;
file3="$1"; shift;

echo "${file1}.png × ${file2}.png -> ${file3}.png"
convert ${file1}.png -colorspace Gray ${file1}.pgm
convert ${file2}.png -colorspace Gray ${file2}.pgm
pnmxarith -multiply ${file1}.pgm ${file2}.pgm > ${file3}.pgm
ls -l ${file1}.pgm ${file2}.pgm ${file3}.pgm
convert ${file3}.pgm -colorspace Gray ${file3}.png
rm ${file1}.pgm ${file2}.pgm ${file3}.pgm
display ${file3}.png
