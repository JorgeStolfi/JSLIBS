#! /bin/bash
# Last edited on 2024-12-24 16:24:22 by stolfi

# Plots a test histogram produced by "pnmxhist"

tmp="/tmp/$$"

hists=( `cd out && ls *.hist | sort | sed -e 's:[.]hist::g'` )

imgs=( `echo ${hists[@]} | tr ' ' '\n' | sed -e 's:-[0-9][0-9]$::g' | sort | uniq` )

colors=( ff0000 008800 0055ff ff8800 008888 8800ff ) # Should be enough...

for coltit in 4:hist 5:cumh ; do
  col="${coltit/:*/}"
  tit="${coltit/*:/}"
  for img in ${imgs[@]} ; do
    chns=( `cd out && ls ${img}-*.hist | sed -e 's:^.*-::g' -e 's:[.]hist::g' | sort | uniq` )
    
    hmin="+inf"
    hmax="-inf"
    # Generate the plot commands: 
    tgfile="${tmp}-${img}.gpl"
    trfile="${tmp}-${img}.range"
    rm -fv ${tgfile}
    printf "plot " > ${tgfile}
    sep='\'
    for c in ${chns[@]}; do 
      # Add the plot command for this channel: 
      hfile="out/${img}-${c}.hist"
      printf "%s\n" "${sep}" >> ${tgfile}
      printf "  \"${hfile}\" using (mid(2,3)):${col} title 'channel ${c}'" >> ${tgfile}
      printf "   with linespoints pt 7 ps 0.75 lc rgb '#${colors[${c}]}'" >> ${tgfile}
      sep=', \'
      
      # Extract and update {hmin,hmax}:
      cat ${hfile} \
        | gawk -v hmin="${hmin}" -v hmax="${hmax}" \
            ' BEGIN { 
                hmin = hmin + 0; hmax = hmax + 0;
              }
              // {
                slo = 0 + $2; shi = 0 + $3;
                if (slo < hmin) { hmin = slo; }
                if (shi > hmax) { hmax = shi; }
              }
              END { printf "%24.16e %24.16e\n", hmin, hmax; }
            ' \
        > ${trfile}
      hrange=( `cat ${trfile}` )
      hmin="${hrange[0]}"
      hmax="${hrange[1]}"
      echo "hmin = ${hmin}  hmax = ${hmax}" 1>&2
    done
    
    # Create the plot file at twice the final scale:
    tpfile="${tmp}-${img}-${c}-${tit}.png"
    GDFONTPATH=ttf
    gnuplot <<EOF
      set term png truecolor rounded font "Arial,28" size 1600,1600
      set output "${tpfile}"
      set title "${img} (${tit})"
      hmin = ${hmin}
      hmax = ${hmax}
      eh = 0.05*(hmax-hmin)
      set xrange [(hmin-eh):(hmax+eh)]
      set yrange [-1: ]
      mid(i,j) = (column(i) + column(j))/2
      load "${tgfile}"
      quit
EOF

    # Shrink the plot file and display it:
    if [[ -s ${tpfile} ]]; then
      pfile="out/${img}-${tit}.png"
      convert ${tpfile} -resize '50%' ${pfile}
      display ${pfile}
    else
      echo "** file ${tpfile} not created" 1>&2; exit 1
    fi
  done
done
