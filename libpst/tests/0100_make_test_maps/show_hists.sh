#! /bin/bash
# Last edited on 2025-02-28 18:44:38 by stolfi

size="$1"; shift

for f in out/*-${size}-Ga{X,Y}-hist.eps ; do
  echo "=== $f ===" 1>&2
  evince $f
done
