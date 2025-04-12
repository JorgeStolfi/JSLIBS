#! /bin/bash
# Last edited on 2025-03-16 19:33:33 by stolfi

ofile="fourob-MF-0512x0422-Z.fni"

tmp="/tmp/$$"

for f in mf-depth mf-confidence ; do 
  cat ${f}.fni \
    | gawk \
        ' /^[ ]*[0-9]+[ ]/ {
            printf "%05d.%05d %s\n", $2, $1, $3
          }
        ' \
    | sort \
    > ${tmp}-${f}.txt
done

join \
    -j1 1 -j2 1 -a 1 -a2 -e '???' \
    -o0,1.2,2.2 \
    ${tmp}-mf-{depth,confidence}.txt \
  > ${tmp}-join.txt

rm ${ofile}
printf "begin float_image_t (format of 2006-03-25)\n" >> ${ofile}
printf "NC = 2\n" >> ${ofile}
printf "NX = 512\n" >> ${ofile}
printf "NY = 422\n" >> ${ofile}

cat ${tmp}-join.txt \
  | gawk \
      ' BEGIN { oy = 0; }
        // { 
          if (NF != 3) { printf "**bug" > "/dev/stderr"; exit(1); }
          xy = $1; dp = $2; wt = $3;
          x = xy; gsub(/^.*[.]/, "", x); x = x + 0;
          y = xy;  gsub(/[.].*$/, "", y); y = y + 0;
          if (dp == "+5.4101562e-01") { 
            dp = "+1.0000000e-05"; wt = "+1.0000000e-05";
          } else {
            wt = 1.0 - wt;
          }
          if (y != oy) { printf "\n"; }
          printf "%5d %5d %s %s\n", x, y, dp, wt;
          oy = y;
        }
      ' \
  >> ${ofile}
      
printf "end float_image_t\n" >> ${ofile}
