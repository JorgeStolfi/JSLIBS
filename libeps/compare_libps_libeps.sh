#! /bin/bash

( cd orig && find ./ -name '*.h' -print ) | sed -e 's:^[.]/::g' > .fnames
( cd orig && find ./ -name '*.c' -print ) | sed -e 's:^[.]/::g' >> .fnames

rm -f .diffs
touch .diffs
for psname in `cat .fnames`; do
  psfile="orig/${psname}"
  epsfile="${psname/pswr/epswr}"
  echo "${psname} ${psfile} ${epsfile}" 1>&2
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
