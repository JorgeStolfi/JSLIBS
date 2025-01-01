#! /bin/bash
# Last edited on 2012-01-27 21:29:59 by stolfilocal

povfile="$1"; shift;

imgfile=${povfile%%.*}.png

if [[ -e /usr/bin/povray ]]; then
  POVRAY="/usr/bin/povray"
elif [[ -e /usr/local/bin/povray ]]; then
  POVRAY="/usr/local/bin/povray"
else
  echo "no povray?" 1>&2; exit 1
fi
POVINC1="/usr/share/povray/include"
POVINC2="povray"
POVTTF="povray/ttf"

${POVRAY} \
    +FN +Q9 +MB1 \
    +W512 +H512 \
    +AM1 +A0.0 +R2 \
    +D +SP32 +EP4 \
    +L${POVINC1} \
    +L${POVINC2} \
    +L${POVTTF} \
    +I${povfile} \
    +O${imgfile} \
  2>&1 | ./povray/povray-output-filter.gawk
display ${imgfile}
