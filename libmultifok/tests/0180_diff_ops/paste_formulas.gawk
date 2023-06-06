#! /usr/bin/gawk -f
# Last edited on 2023-04-28 11:31:59 by stolfi

# Reads one or more formula files as produced by {linear_fit -writeFormula},
# writes them in colums, matching the coeffs by the term name after "#".

# Usage: ${cmd} {file}...

BEGIN { 
  nt = 0; # Number of distinct terms seen.
  split("", termName);  # {termName[0..nt-1]} are the term names.
  split("", termIndex); # {termIndex[na]} is the index {kt} such that {termName[kt]=na}.
  
  nf = 0;
  split("", fname);  # {fname[0..nf-1]} are the names of the files already started processing.
  
  split("", tcoeff); # {tcoeff[kf,kt]} is the coefficient of term {termName[kt]} in file {fname[kf]}.
  
  # Start the term list with the avg and dev errors:
  termName[nt] = "AVG_ERR"; termIndex[termName[nt]] = nt; nt++;
  termName[nt] = "DEV_ERR"; termIndex[termName[nt]] = nt; nt++;

  cur_fname = ""; # Name of current file.
  end_file();
}

/^ *([#]|$)/ {
  next;
}

// {
  fn = FILENAME;
  if (fn != cur_fname) {
    end_file();
    start_file(fn);
  }
}

/^[0-9]+$/ { 
  # Number of terms in file:
  nt_exp = $1 + 0;
  next;
}

/[0-9].*[#]/ {
  if (NF != 4) { data_error("bad coeff line"); }
  it = $1 + 0; cf = $2; cc = $3; na = $4;
  if (it != nt_in_file) { data_error("unexpected term index in line"); }
  if (cc != "#") { data_error("comment char in wrong place"); }
  process_term(na,cf);
  nt_in_file++;
  next;
}

/^ *[-+]?[.0-9]+([e][-+][0-9]+)? *$/ {
  cf = $1;
  if (! has_avg) {
    na = "AVG_ERR"
    has_avg = 1;
  } else if (! has_dev) {
    na = "DEV_ERR"
    has_dev = 1;
  } else {
    data_error("spurious anonymous term line");
  }
  process_term(na,cf);
  next;
}

// { 
  data_error("bad line format");
}

END {
  end_file();
  for (kf = 0; kf < nf; kf++)
    { printf "file %2d = %s\n", kf, fname[kf]; }
  printf "\n"
    
  for (kf = 0; kf < nf; kf++)
    { printf " %16s", ("file " kf); }
  printf " %s", " term name";
  printf "\n"
     
  for (kf = 0; kf < nf; kf++)
    { printf " %16s", "------------"; }
  printf " %s", " ----------------------";
  printf "\n"
 
  for (it = 2; it < nt; it++) {
    na = termName[it];
    print_term(na);
  }
  printf "\n"
  print_term("AVG_ERR");
  print_term("DEV_ERR");
}

function end_file() {
  if (cur_fname != "") {
    printf "file %s - expected %s terms, read %d\n", cur_fname, nt_exp, nt_in_file > "/dev/stderr";
    if (nt_in_file != nt_exp) { data_error("num terms does not match"); }
  }
}

function start_file(fn) {
  fname[nf] = fn;
  nf++;
  nt_in_file = 0; # Number of named term lines read in {cur_fname}.
  nt_exp = -1;    # Number of terms expected in {cur_fname}.
  cur_fname = fn;
  has_avg = 0;    # If 1, the avg error line already appeared.
  has_dev = 0;    # If 1, the error dev line already appeared.
}

function process_term(na,  kt,kf) {
  if (! (na in termIndex)) {
    termName[nt] = na;
    termIndex[na] = nt;
    nt++;
  }
  kt = termIndex[na];
  kf = nf - 1;
  tcoeff[kf,kt] = cf;
}

function print_term(na, kt,cf,kf) {
  kt = termIndex[na];
  for (kf = 0; kf < nf; kf++) {
    if ((kf,kt) in tcoeff) {
      printf " %+16.4f", tcoeff[kf,kt];
    } else {
      printf " %16s", "---";
    }
  }
  printf " %s", na;
  printf "\n";
}

function data_error(msg)
  { printf "%s:%s: ** %s\n", FILENAME, FNR, msg > "/dev/stderr"; 
    printf "  «%s»\n", $0 > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 
          
