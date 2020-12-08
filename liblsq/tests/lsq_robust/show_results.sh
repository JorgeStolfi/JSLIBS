#! /bin/bash
# Last edited on 2014-05-28 12:29:38 by stolfilocal

tag="$1"; shift

if [[ -s out/iter_${tag}_i000.txt ]]; then
  for f in `ls out/iter_${tag}_i*.txt | sort` ; do
    iter="${f##*_i}"; iter="${iter%%.txt}"
    plot_iter.sh iter out/init_${tag}.txt $f &
    sleep 1
  done
fi
plot_iter.sh -1 out/init_${tag}.txt out/term_${tag}.txt

display -title '%f' `ls out/iter_${tag}_i*.png | sort`  out/term_${tag}.png

convert -loop 10 -delay 50 `ls out/iter_${tag}_i*.png | sort` out/iter_${tag}.gif

echo killing `jobs -p`
for p in `jobs -p` ; do /usr/bin/kill -s KILL ${p} ; done
sleep 1
echo did not kill `jobs -p`
