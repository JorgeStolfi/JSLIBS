#! /usr/bin/gawk -f
# Last edited on 2023-02-01 18:40:02 by stolfi

# Creates a file with quadratic term values averaged by sharpness bins.
#
# User must defined (with "-v") the variables 
#
#   {nb}          number of basis elements.
#   {nt}          number of quadratic terms. 
#
# Reads from stdin the "-odata.txt" file produced by {mf_0180_diff_ops.c}. Assumes it has fields
#
#   "P{ki}.{ix}.{iy} {vave} {vdev} {sharp} {zave} {zdev} {coeff[0..nb-1]} {term[0..nt-1]}"
#
# where 
#
#   {ki} is the input image index. 
#   {ix,iy} are column and row of the pixel. 
#   {vave} and {vdev} are the average and deviation of window samples before normalization.
#   {sharp} is the "actual" sharpness of the image at that pixel.
#   {zave} the average in that pixel the of {Z} values relative to {zFoc}.
#   {zdev} is the deviation of {Z} values in that pixel. 
#   {coef[0..nb-1]} are the coefficients of the basis elements in the window samples.
#   {term[0..nt-1]} are the quadratic terms computed from those coefficients.
# 
# Writes to stdout the "-hdata.txt" file with lines 
#
#   "{ih} {sharp_h} {tm_h[0]} ... {tm_h[NT-1]}"
#
# where {ih} is the index and {sharp_h} is the center value of a {sharp} bin 
# and {tm_h[kb]} is the average value of {term[kt]}} for all 
# data lines whose {sharp} that fall in that bin.
#
# A blank line is written if a bin is empty in the middle of the histogram.

BEGIN{ 
  abort = -1;
  if (nb == "") { arg_error("must define {nb}"); } nb += 0;
  if (nt == "") { arg_error("must define {nt}"); } nt += 0;

  split("", term);       # In {term[kt]} is the term of basis element {kt} of one line.

  # In a pass through the input file, we collect {sum_wp_tm[ih][kt]} 
  # which is is the sum of {wp[kp]*term[kp][kt]} for all pixels {kp} 
  # that fall in {sharp} bin {ih}.  We also collect {sum_wp[ih]}.
  split("", sum_wp_tm);    # {sum_wp_tm[ih][kt]} is um of {wp*tm} per histogram bin.
  split("", sum_wp);       # {sum_wp[ih]} is sum of {wp} per histogram bin.
  for (ih = 0; ih < nh; ih++) {
    for (kt = 0; kt < nt; kt++) { sum_wp_tm[ih][kt] = 0; }
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
  nu = 0; # Number of data lines used in histogram.
}

(abort >= 0) { exit(abort); }

/^ *[P][0-9]/ {
  if (NF != 6 + nb + nt) { printf "BUG\n" > "/dev/stderr"; exit(1); }
  id = $1;
  vave = nval($2);
  vdev = nval($3);
  sharp = nval($4); 
  zave = nval($5);
  zdev = nval($6); 
  for (kt = 0; kt < nt; kt++) { term[kt] = nval($(7+nb+kt)); }

  np++;
  
  # wp = sharp/sqrt(sharp*sharp + zdev*zdev); # Any better idea?
  wp = 1.0;

  if (vdev < 0.0) { data_error("bad sample dev"); }
  if (zdev < 0.0) { data_error("bad pixel Z dev"); }
  if (sharp < 0.0) { data_error("bad sample sharpness"); }
  
  accumulate_term_data(vave, vdev, wp, sharp, zave, zdev, nt, term);
  next;
}

// { data_error("invalid line format"); }

END {
  printf "read %d data lines\n", np > "/dev/stderr"
  printf "used %d data lines in histogram\n", nu > "/dev/stderr"
  if (abort >= 0) { exit(abort); }
  write_histogram(nt);
}
  
function accumulate_term_data(vave,vdev,wp,sharp,zave,zdev,nt,term, ih,kt,tm) {
  if ((zdev < 0.3) && (vdev >= 0.02)) {
    ih = floor((sharp - sharpMin)/sharpStep);
    for (kt = 0; kt < nt; kt++) {
      tm = term[kt];
      sum_wp_tm[ih][kt] += wp*tm;
    }
    sum_wp[ih] += wp;
    if(ih > ihmax) { ihmax = ih; }
    if(ih < ihmin) { ihmin = ih; }
    nu++;
  }
}

function write_histogram(nt,  ih,kt,sharp_h,tm_h,nw) {
  nw = 0;
  for (ih = ihmin; ih <= ihmax; ih++) {
    if (sum_wp[ih] > 0) {
      printf "%5d ", ih; 
      sharp_h = sharpMin + (ih + 0.5)*sharpStep;
      printf "%12.6f ", sharp_h; 
      for (kt = 0; kt < nt; kt++) {
        tm_h = sum_wp_tm[ih][kt]/sum_wp[ih]
        printf "  %+14.8f", tm_h;
      }
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

