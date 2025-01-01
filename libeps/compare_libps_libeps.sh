#! /bin/bash

epsdir="/home/stolfi/programs/c/JSLIBS/libeps"
psdir="/home/stolfi/programs/c/JSLIBS-LATER/libps"

( cd ${psdir} && find ./ -name '*.h' -print ) | sed -e 's:^[.]/::g' > .psnames
( cd ${psdir} && find ./ -name '*.c' -print ) | sed -e 's:^[.]/::g' >> .psnames

rm -f .diffs
touch .diffs
for psname in `cat .psnames`; do
  epsname="${psname/pswr/epswr}"
  psfile="${psdir}/${psname}"
  epsfile="${epsdir}/${epsname}"
  echo "${psname} : ${psfile} ${epsfile}" 1>&2
  if [[ -s ${epsfile} ]]; then
    echo "=== ${psfile} ${epsfile} ===" >> .diffs
    echo "" >> .diffs
    cat ${psfile} | sed -f pswr_epswr_rename.sed > .temp
    prdiff -Bb .temp ${epsfile} >> .diffs
  else
    echo "*** No file ${epsfile}"  >> .diffs
  fi
  echo "" >> .diffs
done
