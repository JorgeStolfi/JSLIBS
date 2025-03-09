#! /bin/bash
# Last edited on 2025-02-28 05:57:58 by stolfi

# Plots "${ipref}-cmp.fni" against "${ipref}-ref.fni" in each data channel.

fp="$1"; shift   # "{fnum}-{fname}"
size="$1"; shift # "{nx}x{ny}"
itag="$1"; shift  # 'GS' or 'GE' or 'ZS' or 'ZE'

echo "fp = '${fp}'  size = '${size}'  itag = '${itag}'" 1>&2

ipref="out/${fp}-${size}-${itag}"
tmp="/tmp/$$"

if [[ ( "/${itag}" == "/GS" ) || ( "/${itag}" == "/GE" ) ]]; then
  chans=( 0 1 )
else
  chans=( 0 )
fi

for ch in ${chans[@]}; do 
  # Extract channel ${ch} of each input file:
  cpref="${tmp}-${ch}"
  for sub in "ref" "cmp" ; do
    ifile="${ipref}-${sub}.fni"
    cfile="${cpref}-${sub}.txt"
    echo "creating ${cfile} from channel ${ch} of ${ifile} ..." 1>&2
    cat ${ifile} \
      | gawk -v ch=${ch} \
          ' BEGIN { vmin = +1.0e20; vmax = -1.0e20; }
            /^ *[0-9]+[ ]+[0-9]+/ {
              val = $(ch+3);
              printf "%d:%d %s\n", $1, $2, val; 
              val = val + 0;
              if (val < vmin) { vmin = val; }
              if (val > vmax) { vmax = val; }
            }
            END {
              if (vmin == vmax) { 
                printf "-1:-1 %24.16e\n", vmin - 1.0;
                printf "-2:-2 %24.16e\n", vmax + 1.0;
              }
            }
          ' \
      | sort \
      > ${cfile}
  done

  # Join the two files index by index:
  dfile="${tmp}-ref-cmp.txt"
  echo "joining ${cpref}-ref.txt and ${cpref}-cmp.txt to ${dfile} ..." 1>&2
  join \
      -1 1 -2 1 -o '0,1.2,2.2' -a 1 -a 2 -e '???' \
      ${cpref}-ref.txt ${cpref}-cmp.txt \
    | sort | uniq \
    | tr ':' ' ' \
    | sort -k 3g \
    > ${dfile}

  # Plot one against the other:
  tfile="${tmp}-plot.png"
  echo "plotting ${dfile} to ${tfile} ..." 1>&2
  export GDFONTPATH ttf
  gnuplot <<EOF
    set term png font "arialb,24" size 1800,1600
    set output "${tfile}"
    set xlabel "ref"
    set ylabel "cmp"
    set key top left
    set xzeroaxis
    set yzeroaxis
    plot \
      "${dfile}" using 3:3             title 'ref'   with lines lw 3 lc rgb '#ff8855', \
      "${dfile}" using 3:(2*column(3)) title '2*ref' with lines lw 2 lc rgb '#aaccff', \
      "${dfile}" using 3:4             title 'cmp'   with points pt 7 ps 2.0 lc rgb '#008888'
    quit
EOF
  if [[ -s ${tfile} ]]; then
    pfile="${ipref}-${ch}-ref-cmp-plot.png"
    echo "resizing ${tfile} to ${pfile} ..." 1>&2
    convert ${tfile} -resize '50%' ${pfile}
    display ${pfile}
  else
    echo "** file ${tfile} not generated" 1>&2 ; exit 1
  fi
done
rm -f ${tmp}-*.{txt,png}
