#! /bin/bash
# Last edited on 2013-11-23 01:17:53 by stolfilocal

trial="$1"; shift

if [[ -s out/iter_t${trial}_i000.txt ]]; then
  for f in `ls out/iter_t${trial}_i*.txt | sort` ; do
    plot_iter.sh out/init_t${trial}.txt $f &
    sleep 1
  done
fi
plot_iter.sh out/init_t${trial}.txt out/term_t${trial}.txt

display -title '%f' `ls out/iter_t${trial}_i*.png | sort`  out/term_t${trial}.png

convert -loop 10 -delay 50 `ls out/iter_t${trial}_i*.png | sort` out/iter_t${trial}.gif


echo killing `jobs -p`
for p in `jobs -p` ; do /usr/bin/kill -s KILL ${p} ; done
sleep 1
echo did not kill `jobs -p`
