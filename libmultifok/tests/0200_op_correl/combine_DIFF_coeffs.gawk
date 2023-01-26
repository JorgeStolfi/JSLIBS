#! /usr/bin/gawk -f
# Last edited on 2023-01-25 14:04:36 by stolfi

# Reads the formula output of {lin_fit} for regression with the squared differential 
# operators.  Averages out coeffs of similar terms and prints them in a better order.

# Accepts the formula as written out to {stderr} where each line has "{coef} * X[{it}] # {name}"
# and as written to the formula file, with lines "{it} {coef} # {name}".


BEGIN { 
  split("", coeff) # Indexed by term name, e. g. "DX*DY".
}

/^[ ]*[-+]?[.0-9]+[ ]*[*][ ]*[X]\[[ ]*[0-9]+\][ ]*[#]/ { 
  # A {stderr} tem line.
  gsub(/[ ]*[*] *[X]\[[ 0-9]+\][ ]*[#][ ]*/, " ", $0);
  cf = $1; na = $2;
  printf "%-7s %s\n", na, cf;
  coeff[na] = cf;
  next;
}

/^[ ]*[0-9]+[ ]+[-+]?[.0-9]+[e][-+][0-9]+[ ]*[#]/ { 
  # A formula file term line.
  gsub(/^[ ]*[0-9]+[ ]*/, " ", $0);
  gsub(/[#][ ]*/, " ", $0);
  cf = $1; na = $2;
  printf "%-7s %s\n", na, cf;
  coeff[na] = cf;
  next;
}

//{ next }

END {
  pt1("F*F");
  pt2("DX*DX", "DY*DY");
  pt2("DXX*DXX", "DYY*DYY");
  pt1("DXY*DXY");
  pt2("S3*S3", "C3*C3");
  pt1("Q*Q");
  pt1("1");
}

function pt1(na) {
  if (na in coeff) { 
    printf "%-8s %+12.4f\n", na, coeff[na];
  }
}

function pt2(na1,na2, cav) {
  if ((na1 in coeff) || (na2 in coeff)) {
    cav = (coeff[na1] + coeff[na2])/2;
    printf "%-8s %+12.4f\n", na1, cav;
    printf "%-8s %+12.4f\n", na2, cav;
  }
}
