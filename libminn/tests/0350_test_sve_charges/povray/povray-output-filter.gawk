#! /usr/bin/gawk -f
# Last edited on 2011-04-22 22:30:37 by stolfi

BEGIN { 
  state = 0;  
    # If 0, just copying stream.
    # If 1, deleting irrelevant output.
    # If 2, saving contex of error message in {ectx[0..nctx-1]}.
    # If 3, collecting continuation lines of error message text {emess}. 
}
 
/^[ ]/ { 
  # Possible continuation line:
  if (state == 3) { emess = (emess $0); next; }
}

// { 
  # Not a continuation line. If gathering error message, dump it:
  if (state == 3) { 
    if ($0 ~ /^[ ]/) { prog_error(("boh?")); }
    dump_povray_message()
    state = 0;
  }
}

/^Persistence of Vision[(]tm[)]/ { state = 1; next; }
/^Output Options/ { state = 0; }
/^Total Scene Processing Times/ { state = 1; next; }

/^File[:] *[^ ]*[ ]+Line[:] */ {
  # Preamble of context lines and error message:
  efile=$2; eline=$4;
  emess="??";
  split("", ectx); nctx = 0;
  state = 2; 
  next;
}

/^File Context [(][0-9]+[ ]+lines[)][:]/ {
  # Header of context lines, ignore:
  if (state != 2) { prog_error(("seq bug")); }
  next;
}

/^Parse Error/ { 
  # End of context lines, error message (may continue):
  if (state != 2) { prog_error(("seq bug")); }
  emess = $0; gsub(/Parse Error: /, "", emess); 
  state = 3; 
  next;
}

/^[!][!]/ { 
  # User debugging line:
  print > "/dev/stderr";
  next;
}

(state == 0) {
  # Generic random line:
  print > "/dev/stderr"; 
  next;
}

(state == 1) { 
  # Undesired noise line:
  next;
}

(state == 2) {
  # Error context line, save it:
  ectx[nctx] = ("  " $0); nctx++;
  next;
}

(state == 3) {
  # Should have been handled before:
  prog_error(("boh?")); 
}

END {
  if (state == 3) { 
    dump_povray_message()
    state = 0;
  }
  printf "\n" > "/dev/stderr"; 
  printf "--- povray finished ------------------------------------------------\n" > "/dev/stderr"; 
}

function dump_povray_message(   i) {
  printf "%s:%s: ** %s\n", efile, eline, emess > "/dev/stderr";
  for (i = 0; i < nctx; i++) { printf "%s\n", ectx[i] > "/dev/stderr"; }
  nctx = 0;
}

function prog_error(msg) {
  printf "!! povray-output-filter: bug at output line %d - %s\n", FNR, msg > "/dev/stderr";
  print > "/dev/stderr";
  next;
}
