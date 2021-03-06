#! /bin/bash
# Last edited on 2017-03-13 22:35:18 by stolfilocal

# Plots the curve generated by {test_sve_catenary}

file="$1"; shift;

tmp=/tmp/$$

if [[ -r ${file} ]] ; then
  fni_plot.sh < ${file} > ${tmp}.eps
  display -flatten ${tmp}.eps
  rm -f ${tmp}.eps
else
  echo "file ${file} missing" 1>&2
fi
