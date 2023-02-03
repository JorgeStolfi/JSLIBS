#! /usr/bin/gawk -f
# Last edited on 2023-02-01 19:02:02 by stolfi

BEGIN {
  np = 0;          # Number of records read.

  zdev_max = 1.0;  # Reject pixels with {zdev} above  this.
  nr_zdev = 0;     # Number of records rejected for high {zdev}.
  
  sharp_min = 0.4; # Reject pixels with {sharp} below this.
  nr_sharp = 0;    # Number of records rejected for low {sharp}.

  nu = 0;          # Number of records kept.
}

// {
  # Remove commants:
  gsub(/[ ]*[#].*$/, "", $0);
}

/^ *$/ { 
  # Blank line, igore.
  next;
}

/^ *P[0-9]/{ 
  if (NF < 6) {  printf "%s:%d: ** wrong number of fields = %d\n", FILENAME, FNR, NF > "/dev/stderr"; exit(1); }
  ID = $1; vavg = $2+0; vdev = $3+0; sharp = $4+0; zave = $5+0; zdev = $6+0; 
  np++;
  # Filter out bad pixel data:
  if (sharp < sharp_min) { nr_sharp++; next; }
  if (zdev > zdev_max) { nr_zdev++; next; }
  print;
  nu++;
  next
  next
}


// { printf "%s:%d: ** bad pixel data line format\n", FILENAME, FNR > "/dev/stderr"; exit(1); }


END {
  printf "filter_pixel_data.gawk:\n" > "/dev/stderr";
  printf "  read %d pixel data lines\n", np > "/dev/stderr";
  printf "  rejected %d lines for low {sharp}\n", nr_sharp > "/dev/stderr";
  printf "  rejected %d lines for high {zdev}\n", nr_zdev > "/dev/stderr";
  printf "  ueed %d lines\n", nu > "/dev/stderr";
}
