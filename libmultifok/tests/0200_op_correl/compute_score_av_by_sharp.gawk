#! /usr/bin/gawk -f
# Last edited on 2023-01-23 11:39:48 by stolfi

# Used by {plot_regression.sh}.
# Reads the join of "-cdata.txt" and "-regr.txt" files.
# Outputs a file with a selected {score} column averaged 
# by bins of actual {sharp}.

# User must define (with "-v") the variables 
# {col}      index of score colum to process.
# {sqrSharp} If 1, assumes that {sharp} and {score} are squared in the regression file. 
# {nh}       number of {sharp} bins.

BEGIN { 
  if (sqrSharp == "") { arg_error("must define {sqrSharp}"); } sqrSharp += 0;
  if (nh == "") { arg_error("must define {nh}"); } nh += 0;
  if (col == "") { arg_error("must define {col}"); } col += 0;
  split("", sum_wp_sc); 
  split("", sum_wp); 
  ihmin = +9999;
  ihmax = -9999;
  nr = 0 # Records read.
}
/^ *P/ { 
  id = $1; wp = nval($2); sharp_c = nval($3); score_c = nval($4); sharp_r = nval($5); score_r = nval($6);
  nr++;
  if (wp <= 0.0) { data_error("bad weight"); }
  if (sqrSharp)
    { # Un-square the {sharp} and {score} from the {regr} file for compatibility:
      sharp_r = sqrt(sharp_r); 
      score_r = (score_r < 0.0 ? -sqrt(-score_r) : sqrt(score_r));
    }
  if (fabs(sharp_c - sharp_r) > 1.0e-5) 
    { data_error(("join error " sharp_c " " sharp_r)); }
  sharp = sharp_r;  # May be squared.
  if (col == 4) {
    score = score_c;
  } else if (col == 6) {
    score = score_r;
  } else {
    arg_error(("bad {col} " col)); 
  }
  
  ih = int(nh*sharp);
  if (! (ih in sum_wp)) { 
    sum_wp[ih] = 0;
    sum_wp_sc[ih] = 0;
    if (ih < ihmin) { ihmin = ih; }
    if (ih > ihmax) { ihmax = ih; }
  }
  sum_wp_sc[ih] += wp*score; 
  sum_wp[ih] += wp;
}

END {
  if (ihmin > ihmax) { data_error("no data?"); }
  for (ih = ihmin; ih <= ihmax; ih++) { 
    if (ih in sum_wp) { 
      score_av = sum_wp_sc[ih]/sum_wp[ih]; 
      sharp_av = (ih + 0.5)/nh;
      printf "%5d %16.12f %16.12f\n", ih, sharp_av, score_av;
    }
  }
}

function fabs(x) {
  return (x < 0 ? -x : x);
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


