#! /usr/bin/gawk -f
# Last edited on 2023-01-25 14:12:08 by stolfi

# Reads "-odata.txt" file produced by {mf_0180_diff_ops.c}
# Writes a file with lines 
#   "{ih} {sharp_h} {wr_h} {cfq_h[0]} ... {cfq_h[NB-1]}" for plotting.
# where {ih} is the index and {sharp_h} is the center value of a {sharp} bin, 
# {wr_h} is a bin weight for regression, 
# and {cfq_h[kb]} is the average value of {coeff[kb]/vdev} squared for 
# data lines that fall in that bin.
#
# If {unitTerm} parameter is true (1), an extra " 1.0" is appended after the last coefficient.
# 
# Blank lines are written between images.

BEGIN{ 
  abort = -1;
  if (nb == "") { arg_error("must define {nb}"); } nb += 0;
  if (unitTerm == "") { arg_error("must define {unitTerm}"); } unitTerm += 0;

  split("", coeff);       # In {coeff[kb]} is the coeff of basis element {kb} of one line.

  # In a pass through the input file, we collect {sum_wp_cfq[ih][kb]} 
  # which is is the sum of {wp[kp]*coeff[kp][kb]} for all pixels {kp} 
  # that fall in {sharp} bin {ih}.  We also collect {sum_wp[ih]}.
  split("", sum_wp_cfq);   # Sum of {wp*cfq} per histogram bin.
  split("", sum_wp);       # Sum of {wp} per histogram bin.
  for (ih = 0; ih < nh; ih++) {
    for (kb = 0; kb < nb; kb++) { sum_wp_cfq[ih][kb] = 0; }
    sum_wp[ih] = 0;
  }
 
  # Parameters of histogram:
  sharpMin = 0.0;
  sharpMax = 1.0;
  sharpStep = 0.02;
  nh = ceil((sharpMax - sharpMin)/sharpStep) + 1; # Max number of hist bins.
  ihmin = +9999; # Index of first non-empty bin. 
  ihmax = -9999; # Index of last non-empty bin. 
  np = 0; # Number of data lines read.
  nu = 0; # Number of pixels used in histogram.
}

(abort >= 0) { exit(abort); }

/^ *[P][0-9]/ {
  id = $1;
  vave = nval($2);
  vdev = nval($3);
  wp = nval($4); 
  sharp = nval($5); 
  zave = nval($6);
  zdev = nval($7); score_r = nval($6);
  for (kb = 0; kb < nb; kb++) { coeff[kb] = nval($(8+kb)); }

  np++;

  if (wp <= 0.0) { data_error("bad pixel weight"); }
  if (vdev < 0.0) { data_error("bad sample dev"); }
  if (sharp < 0.0) { data_error("bad sample sharpness"); }
  
  accumulate_coeff_data(vave, vdev, wp, sharp, zave, zdev, nb, coeff);
  next;
}

// { data_error("invalid line format"); }

END {
  if (abort >= 0) { exit(abort); }
  write_histogram(nb);
}
  
function accumulate_coeff_data(vave,vdev,wp,sharp,zave,zdev,nb,coeff, ih,kb,cfq) {
  if ((zdev < 0.3) && (vdev >= 0.02)) {
    ih = floor((sharp - sharpMin)/sharpStep);
    for (kb = 0; kb < nb; kb++) {
      cfq = coeff[kb]*coeff[kb];
      sum_wp_cfq[ih][kb] += wp*cfq;
    }
    sum_wp[ih] += wp;
    if(ih > ihmax) { ihmax = ih; }
    if(ih < ihmin) { ihmin = ih; }
    nu++;
  }
}

function write_histogram(nb,  ih,kp,sharp_h,wr_h,cfq_h,nw) {
  nw = 0;
  for (ih = ihmin; ih <= ihmax; ih++) {
    if (sum_wp[ih] > 0) {
      printf "%5d ", kb; 
      sharp_h = sharpMin + (ih + 0.5)*sharpStep;
      printf "%+12.6f ", sharp_h; 
      wr_h = sharp_h;
      printf "%8.4f ", wr_h; 
      for (kb = 0; kb < nb; kb++) {
        cfq_h = sum_wp_cfq[ih][kb]/sum_wp[ih]
        printf "  %12.6f", cfq_h;
      }
      if (unitTerm) { printf "  1.0"; }
      nw++
    }
    printf "\n"; # If no data in bin, print a blank line to break the plot.
  }
  printf "\n";
  printf "Read %d lines, used %d, %d discarded", np, nu, np-nu > "/dev/stderr"; 
  printf " %d plot points\n", nw > "/dev/stderr";
}
 
function fabs(x) {
  return (x < 0 ? -x : x);
}

function floor(x)
  { x = (x >= 0? int(x) : -int(-x));
    return x;
  }
 
function ceil(x)
  { x = (x <= 0? int(x) : -int(-x));
    return x;
  }

function nval(x) {
  if (x == "??") { data_error("unpaired line"); }
  return x + 0;
}
         
function data_error(msg)
  { printf "%s:%s: ** %s\n", FILENAME, FNR, msg > "/dev/stderr"; 
    printf "  «%s»\n", $0 > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 
          
function arg_error(msg)
  { printf "** %s\n", msg > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 

