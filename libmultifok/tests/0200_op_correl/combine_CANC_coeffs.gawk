#! /usr/bin/gawk -f
# Last edited on 2023-01-23 04:23:52 by stolfi

# Reads the formula output of {lin_fit} for regression with all terms of order 0 and 2
# from the canonical basis.  Averages out coeffs of terms that should have same weight
# because of rotation and flip symmetries, and prints them in a better order.

# Accepts the formula as written out to {stderr} where each line has "{coef} * X[{it}] # {name}"
# and as written to the formula file, with lines "{it} {coef} # {name}".

BEGIN { 
  split("", coeff) # Indexed by term name, e. g. "Pmo*Pop".
}

/^[ ]*[-+]?[.0-9]+[ ]*[*][ ]*[X]\[[ ]*[0-9]+\][ ]*[#]/ { 
  # A {stderr} tem line.
  gsub(/[ ]*[*] *[X]\[[ 0-9]+\][ ]*[#][ ]*/, " ", $0);
  cf = $1; term = $2;
  printf "%-7s %s\n", term, cf "/dev/stderr";
  coeff[term] = cf;
  next;
}

/^[ ]*[0-9]+[ ]+[-+]?[.0-9]+[e][-+][0-9]+[ ]*[#]/ { 
  # A formula file term line.
  gsub(/^[ ]*[0-9]+[ ]*/, " ", $0);
  gsub(/[#][ ]*/, " ", $0);
  cf = $1; term = $2;
  printf "%-7s %s\n", term, cf > "/dev/stderr";
  coeff[term] = cf;
  next;
}

//{ next }

END {
  split("", sum_coeff); # Indexed by reduced term, sum of coeffs of all related terms.
  split("", num_coeff); # Indexed by reduced term, count of all related terms.
  split("Poo Pom Pop Pmo Ppo Pmm Ppp Pmp Ppm", els); # Basis element names.
  
  printf "generating terms...\n" > "/dev/stderr";
  split("", all_terms); # All terms that may appear in the data, from 0 to {nt-1}.
  split("", term_index); # Maps term name to term index.
  split("", red_term);  # Maps each term name to its reduced name.
  nt = make_all_terms(els, all_terms, red_term); 
  for (kt = 0; kt < nt; kt++) { 
    term = all_terms[kt]
    printf "  %-8s = %s\n", term, red_term[term] > "/dev/stderr";
    term_index[term] = kt;
  }
 
  printf "combining coeffs of equivalent terms...\n" > "/dev/stderr";
  for (term in coeff) {
    if (! (term in term_index)) { prog_error(("duh? " term)); }
    rterm = red_term[term];
    sum_coeff[rterm] += coeff[term];
    num_coeff[rterm] ++;
  }
 
  printf "printing equalized coeffs...\n" > "/dev/stderr";
  for (kt = 0; kt < nt; kt++) {
    term = all_terms[kt];
    rterm = red_term[term];
    seq = sum_coeff[rterm];
    neq = num_coeff[rterm];
    if (neq > 0) {
      cf = seq/neq
      printf "  %-8s = %-8s %+12.4f -> %+14.8f / %3d %+12.4f\n", term, red_term[term], coeff[term], seq, neq, cf  > "/dev/stderr";
      printf "%-8s %+12.4f\n", term, cf;
    }
  }
}

function make_all_terms(els,at,rna, k1,k2,el1,el2,kt,na) {
  kt = 0;
  for (k1 in els) {
    el1 = els[k1];
    for (k2 in els) {
      el2 = els[k2];
      if (el1 <= el2) {
        na = (el1 "*" el2)
        at[kt] = na;
        rna[na] = reduced_term(el1, el2);
        kt++;
      }
    }
  }
  at[kt] = "1"; kt++;
  rna["1"] = "1";
  return kt;
}

function reduced_term(el1,el2,  ks,kv,kh,eterm,rterm) {
  if (el1 !~ /^P[mop][mop]$/) { prog_error(("bad basis elem name {el1} = " el1)); }
  if (el2 !~ /^P[mop][mop]$/) { prog_error(("bad basis elem name {el2} = " el2)); }
  rterm = "";
  for (ks = 0; ks < 2; ks++) {
    for (kv = 0; kv < 2; kv++) {
      for (kh = 0; kh < 2; kh++) {
        if (el1 <= el2) {
          eterm = (el1 "*" el2);
        } else {
          eterm = (el2 "*" el1);
        }
        if ((rterm == "") || (eterm < rterm)) { rterm = eterm; }
        el1 = flip_hor(el1);
        el2 = flip_hor(el2);
      }
      el1 = flip_ver(el1);
      el2 = flip_ver(el2);
    }
    el1 = swap(el1);
    el2 = swap(el2);
  }
  return rterm;
}

function flip_hor(el, x,y) {
  x = substr(el,2,1);
  y = substr(el,3,1);
  return ("P" (x == "m" ? "p" : (x == "p" ? "m" : x)) y);
}

function flip_ver(el, x,y) {
  x = substr(el,2,1);
  y = substr(el,3,1);
  return ("P" x (y == "m" ? "p" : (y == "p" ? "m" : y)));
}

function swap(el, x,y) {
  x = substr(el,2,1);
  y = substr(el,3,1);
  return ("P" y x);
}
         
function data_error(msg) { 
  printf "%s:%s: ** %s\n", FILENAME, FNR, msg > "/dev/stderr"; 
  printf "  «%s»\n", $0 > "/dev/stderr"; 
  abort = 1;
  exit(abort);
} 
          
function arg_error(msg) { 
  printf "** %s\n", msg > "/dev/stderr"; 
  abort = 1;
  exit(abort);
} 
          
function prog_error(msg) { 
  printf "** prog error: %s\n", msg > "/dev/stderr"; 
  abort = 1;
  exit(abort);
} 
