#! /usr/bin/gawk -f
# Last edited on 2023-04-28 11:34:12 by stolfi

# Creates a file with quadratic term values averaged by bins of sampling radius squared.
#
# User must defined (with "-v") the variables 
#
#   {nb}          number of basis elements.
#   {nt}          number of quadratic terms. 
#
# Reads from stdin the "-odata.txt" file produced by {mf_0180_diff_ops.c}. Assumes it has fields
#
#   "P{ki}.{ix}.{iy} {vave} {vgrd} {vdev} {sharp} {zrav} {zdev} {coeff[0..nb-1]} {term[0..nt-1]}"
#
# where 
#
#   {ki} is the input image index. 
#   {ix,iy} are column and row of the pixel. 
#   {vave} is the average of window samples before normalization.
#   {vgrd} is the gradient modulus in window before normalization.
#   {vdev} is the rms of window samples after removing average and gradient.
#   {sharp} is the "actual" sharpness of the image at that pixel.
#   {zrav} the average in that pixel the of {Z} values relative to {zFoc}.
#   {zdev} is the deviation of {Z} values in that pixel. 
#   {coef[0..nb-1]} are the coefficients of the basis elements in the window samples.
#   {term[0..nt-1]} are the quadratic terms computed from those coefficients.
# 
# Writes to stdout the "-hdata.txt" file with lines 
#
#   "{i} {hrad_h[i]} {tm_h[i][0]} ... {tm_h[i][NT-1]}"
#
# where {i} is the index of a bin, and {hrad_h[i]} is the central value of 
# the sampling radius in the bin, and {tm_h[i][kb]} is the average value of {term[kt]}} for all 
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
  # that fall in bin {ih}.  We also collect {sum_wp[ih]}.
  split("", sum_wp_tm);    # {sum_wp_tm[ih][kt]} is um of {wp*tm} per histogram bin.
  split("", sum_wp);       # {sum_wp[ih]} is sum of {wp} per histogram bin.
  for (ih = 0; ih < nh; ih++) {
    for (kt = 0; kt < nt; kt++) { sum_wp_tm[ih][kt] = 0; }
    sum_wp[ih] = 0;
  }
 
  # Parameters of histogram:
  hradMin = 0.0;
  hradMax = 30.0;
  hradStep = 0.25;
  nh = ceil((hradMax - hradMin)/hradStep) + 1; # Max number of hist bins.
  ihmin = +9999; # Index of first non-empty bin. 
  ihmax = -9999; # Index of last non-empty bin. 
  np = 0; # Number of data lines read.
}

(abort >= 0) { exit(abort); }

/^ *[P][0-9]/ {
  if (NF != 7 + nb + nt) { printf "BUG\n" > "/dev/stderr"; exit(1); }
  id = $1;
  vave = nval($2);
  vgrd = nval($3);
  vdev = nval($4);
  sharp = nval($5); 
  zrav = nval($6);
  zdev = nval($7); 
  
  if (sharp <= 0.0) { data_error("bad sharpness"); }
  hrad = 1.0/sharp;
  for (kt = 0; kt < nt; kt++) { term[kt] = nval($(8+nb+kt)); }

  np++;
  
  # wp = sharp/sqrt(sharp*sharp + zdev*zdev); # Any better idea?
  wp = 1.0;

  if (vdev < 0.0) { data_error("bad sample dev"); }
  if (zdev < 0.0) { data_error("bad pixel Z dev"); }
  
  accumulate_term_data(vave, vdev, wp, hrad, zrav, zdev, nt, term);
  next;
}

// { data_error("invalid line format"); }

END {
  printf "average_terms_by_sharp.gawk:\n" > "/dev/stderr";
  printf "  read %d data lines\n", np > "/dev/stderr"
  printf "  bin index range %d..%d (%d bins)\n", ihmin, ihmax, ihmax-ihmin+1 > "/dev/stderr"
  if (abort >= 0) { exit(abort); }
  write_histogram(nt);
}
  
function accumulate_term_data(vave,vdev,wp,hrad,zrav,zdev,nt,term, ih,kt,tm) {
  
  ih = floor((hrad - hradMin)/hradStep);
  for (kt = 0; kt < nt; kt++) {
    tm = term[kt];
    sum_wp_tm[ih][kt] += wp*tm;
  }
  sum_wp[ih] += wp;
  if(ih > ihmax) { ihmax = ih; }
  if(ih < ihmin) { ihmin = ih; }
}

function write_histogram(nt,  ih,kt,hrad_h,tm_h,nw) {
  nw = 0;
  for (ih = ihmin; ih <= ihmax; ih++) {
    if (sum_wp[ih] > 0) {
      printf "%5d ", ih; 
      hrad_h = hradMin + (ih + 0.5)*hradStep;
      printf "%12.6f ", hrad_h; 
      for (kt = 0; kt < nt; kt++) {
        tm_h = sum_wp_tm[ih][kt]/sum_wp[ih]
        printf "  %+14.8f", tm_h;
      }
      nw++
    }
    printf "\n"; # If no data in bin, print a blank line to break the plot.
  }
  printf "\n";
  printf "  %d non-empty bins written to histogram file\n", nw > "/dev/stderr";
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

