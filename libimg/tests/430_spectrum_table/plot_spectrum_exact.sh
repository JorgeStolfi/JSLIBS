#! /bin/bash
# Last edited on 2023-03-19 14:22:11 by stolfi

# Plots the "exact" spectrum table produced by {test_spectrum_table}.
# Usage: "plot_spectrum_exact {NAME}"
# Reads file "{NAME}.txt", writes "{NAME}-n.eps", "{NAME}-p.eps"

name="$1"; shift;

# Convert the exact spectrum "${name}.txt" to a cumulative histogram "${tmp}.txt".
# For each distinct absolute frequency {FREQ} in the input, the output contains
# two lines
#   "{FREQ} -1 {NTERMS_LSS} {POWER_LSS}" 
#   "{FREQ} 00 {NTERMS_EQL} {POWER_EQL}" 
#   "{FREQ} +1 {NTERMS_LEQ} {POWER_LEQ}" 
# The fields {NTERMS_XXX} and {POWER_XXX} are the sums of {NTERMS}
# and {POWER} for the following sets of exact spectrum entries:
#   "{XXX}" = "{LSS}":  absolute frequency *less than* {FREQ};
#   "{XXX}" = "{EQL}":  absolute frequency *equal to* {FREQ};
#   "{XXX}" = "{LEQ}":  absolute frequency *less than or equal to* {FREQ};
#
# Makes sure that the first line is 
#   "0.0000 -1 0.000 0.000"
# and the last two lines are
#   "{sqrt(0.5)} +1 {TOT_NTERMS} {TOT_POWER}"
#   "{sqrt(0.5)} -1 0.000 0.000"

tmp="/tmp/$$"

cat ${name}.txt \
  | gawk \
      ' BEGIN { 
          ofr = 0; # {FREQ} of previous record.
          tn_lss = 0; # Total {NTERMS} for terms with {FREQ < ofr}
          tp_lss = 0; # Total {POWER} for terms with {FREQ < ofr}
          start_new_freq(fr);
        } 
        /[0-9]/ { 
          f2n = $1; f2d = $2; fr = $3; n = $4; p = $5; 
          if (fr < ofr) { printf "freq out of order\n" > "/dev/stderr"; exit(1); }
          if (fr != ofr) { 
            finish_prev_freq();
            start_new_freq(fr);
          }
          # Accumulate this term into {tn_eql,tp_eql}:
          tn_eql += n; tp_eql += p;
        }
        function finish_prev_freq()
        { # Output the {EQL} and {LEQ} lines for frequency {ofr}:
          if ((tn_eql != 0) || (tp_eql != 0)) 
            { printf "%24.16e  00 %24.16e %24.16e\n", ofr, tn_eql, tp_eql; }
          printf "%24.16e  +1 %24.16e %24.16e\n", ofr, tn_lss + tn_eql, tp_lss + tp_eql;
          # Update {tn_lss,tp_lss} expecting the next freq:
          tn_lss += tn_eql; tp_lss += tp_eql;
        }
        function start_new_freq(fr)
        { # Reset {tn_eql,tp_eql} and update {ofr} for the new freq {fr}: 
          tn_eql = 0; tp_eql = 0;
          ofr = fr;
          # Output the {LSS} line for the new frequency {ofr}:
          printf "%24.16e  -1 %24.16e %24.16e\n", ofr, tn_lss, tp_lss;
        }
        END {
          # Force the plot to end with {FREQ == sqrt(0.5)}: 
          fr = sqrt(0.5);
          if (fr != ofr) { 
            finish_prev_freq();
            start_new_freq(fr);
          }
          finish_prev_freq();
          # Dummy line to force the bar plot to come back to the axis:
          tn_lss = 0; tp_lss = 0;
          start_new_freq(fr);
        }
      ' \
  > ${tmp}.txt

# Get the total number of terms:
tot_terms=`cat ${tmp}.txt | tail -2 | head -1 | gawk '//{printf "%d\n", $3+0;}'`
echo "tot_terms = ${tot_terms}"
        
# Plot the cumulative histograms:
gnuplot <<EOF

set term postscript eps color lw 2 "TimesRoman" 24
set size 1,1
set yrange [-0.00000001:]
set xrange [-0.023:+0.730]
set nokey

# Area of Maltese cross for abs freq {f} (waves/pixel):
across(f) = 2*sqrt(f**2 - 0.25)
# Area of circular sectors betwen arms of Maltese cross:
asects(f) = f**2*(pi - 4*acos(1/(2*f)))
# Total area:
area(f) = (f < 0.5 ? pi*f**2 : (f < sqrt(2) ? across(f)+asects(f) : 1.0))

set title "exact spectrum - terms"
set output "${name}-in.eps"
plot "< ./pick_impulses.sh < ${tmp}.txt" using 1:(column(3)+0.00001) with impulses

set title "exact spectrum - power"
set output "${name}-ip.eps"
plot "< ./pick_impulses.sh < ${tmp}.txt" using 1:(column(4)+0.00001) with impulses

set title "exact spectrum - cumulative terms"
set output "${name}-cn.eps"
plot \
  ((${tot_terms})*area(x)) with lines lt 2, \
  "< ./pick_steps.sh < ${tmp}.txt" using 1:(column(3)+0.00001) with lines lt 1

set title "exact spectrum - cumulative power"
set output "${name}-cp.eps"
plot "< ./pick_steps.sh < ${tmp}.txt"  using 1:(column(4)+0.00001) with lines

set title "exact spectrum - cum power/terms"
set output "${name}-cr.eps"
plot "< ./pick_steps.sh < ${tmp}.txt"  using 1:(column(4)/column(3)) with lines

EOF
 
rm ${tmp}.txt

