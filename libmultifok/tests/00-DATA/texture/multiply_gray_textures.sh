#! /bin/bash
# Last edited on 2025-04-19 08:10:14 by stolfi

# Multiplies grayscale images {file1} and {file2}, cropped
# to {NX} by {NY}, and writes the result as {file3}.

file1="$1"; shift;
file2="$1"; shift;
NX="$1"; shift
NY="$1"; shift
file3="$1"; shift;

echo "${file1}.png × ${file2}.png -> ${file3}.png (${NX} × ${NY})"

tmp="/tmp/$$"
pgm1="${tmp}-1"; shift;
pgm2="${tmp}-1"; shift;
pgm3="${tmp}-3.pgm"; shift;

convert ${file1} -colorspace Gray -crop '${NX}x${NY}' ${pgm1}
convert ${file2} -colorspace Gray -crop '${NX}x${NY}' ${pgm2}
pnmxarith -multiply ${pgm1} ${pgm2} > ${pgm3}
identify ${pgm1} ${pgm2} ${pgm3}
convert ${pgm3} -colorspace Gray ${file3}
rm -fv ${pgm1} ${pgm2} ${pgm3}
display ${file3}
