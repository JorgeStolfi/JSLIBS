#! /bin/bash
# Last edited on 2023-03-29 19:38:15 by stolfi

(cd out && ls -tr *-euler-eerr.plot ) \
  | sed -e 's/-euler-eerr.plot//g' \
  > out/.problems

(cd out && ls -tr quadr-*-eerr.plot ) \
  | sed -e 's/^quadr-//g' -e 's/-eerr.plot//g' \
  > out/.methods

for prbl in `cat out/.problems` ; do
  for mthd in `cat out/.methods` ; do
    show_error_plot.sh ${prbl} ${mthd} 15
  done
done

