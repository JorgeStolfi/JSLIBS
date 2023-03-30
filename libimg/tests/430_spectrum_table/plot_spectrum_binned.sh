#! /bin/bash
# Last edited on 2023-03-19 14:22:28 by stolfi

# Plots the "binned" spectrum table produced by {test_spectrum_table}.
# Usage: "plot_spectrum_binned {NAME}"
# Reads file "{NAME}.txt", writes "{NAME}-n.eps", "{NAME}-p.eps"

name="$1"; shift;
 
# Convert the binned spectrum "${name}.txt" to a histogram plot "${tmp}.txt".
# Each line has "{FREQ} {TOT_NTERMS} {TOT_POWER}" 

tmp="/tmp/$$"
cat ${name}.txt \
  | gawk \
      ' BEGIN { ofr = 0; on = 0; op = 0;
          printf "%24.16e  %24.16e %24.16e\n", ofr, on, op; 
        } 
        /[0-9]/ { 
          fmin = $1; fmax = $2; fr = $3; n = $4; p = $5; 
          if (fmin < ofr) { printf "freq out of order\n" > "/dev/stderr"; exit(1); }
          if (ofr < fmin) { 
            on = 0; op = 0;
            printf "%24.16e  %24.16e %24.16e\n", ofr, on, op;
          }
          ofr = fmin; on = n; op = p;
          printf "%24.16e  %24.16e %24.16e\n", ofr, on, op; 
          ofr = fmax;
          printf "%24.16e  %24.16e %24.16e\n", ofr, on, op; 
        }
        END {
          fmin = sqrt(0.5); on = 0; op = 0;
          if (ofr < fmin) { 
            printf "%24.16e  %24.16e %24.16e\n", ofr, on, op;
          }
          ofr = fmin; 
          printf "%24.16e  %24.16e %24.16e\n", ofr, on, op;
        }
      ' \
  > ${tmp}.txt
       
# Plot the histograms:
gnuplot <<EOF

set term postscript eps color lw 2 "TimesRoman" 24
set size 1,2
set yrange [-0.00000001:]
set xrange [-0.023:+0.730]
set nokey

set title "binned spectrum - bins"
set output "${name}-b.eps"
plot "${tmp}.txt" using 1:(1 + (int(column(0)) % 4 - 1.5)**2) with lines

set title "binned spectrum - terms"
set output "${name}-n.eps"
plot "${tmp}.txt" using 1:2 with lines

set title "binned spectrum - power"
set output "${name}-p.eps"
plot "${tmp}.txt"  using 1:3 with lines

set yrange [+0.0001:]
set logscale y
set title "binned spectrum - power/terms"
set output "${name}-r.eps"
plot "${tmp}.txt"  using 1:(column(3)/column(2)) with lines

EOF
 
