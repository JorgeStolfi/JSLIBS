#! /usr/bin/gawk -f
# Last edited on 2002-03-16 03:35:09 by stolfi

# Converts a list of "xv" mouse click data to a seed file, as expected
# by "pnmift".
# 
# If the "scale" variable is set, multiplies it into the coordinates.
#
# All lines beginning with "#" lines are preserved. Lines beginning
# with "##" define the label to use for the following clicks.

/^[#][#]/ {
  label = $2;
}

/^[#]/ {
  print; 
  next;
}

/^[ ]*$/ {
  print "# ";
  next;
}

/^[ ]*[(][ ]*[0-9]+[ ]+[0-9]+[ ]*[)][ ]*[=][ ]*[0-9]+[ ]*$/ {
  # Line was already cleaned - preserve coordinates, but change label:
  gsub(/[()=]/, " ", $0);
  col = $1;
  row = $2;
  lab = $3;
  printf "( %4d %4d ) = %5d\n", col, row, label; 
  next;
}

/^[ ]*[0-9]+[,][ ]*[0-9]+[ ]*[=]/ {
  # xv mouse click - reformat:
  col = $1; gsub(/[,]*$/, "", col);
  row = $2; gsub(/[=]*$/, "", row);
  lab = label;
  if (scale != "")
    { col = int(scale * col + 0.5); 
      row = int(scale * row + 0.5); 
    }
  printf "( %4d %4d ) = %5d\n", col, row, lab; 
  next;
}
  
/./{
  printf "line %d: bad click format\n", FNR > "/dev/stderr";
  exit 1;
}
