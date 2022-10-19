#! /usr/bin/gawk -f
# Last edited on 2021-12-19 07:27:41 by stolfi

# Reads from {stdin} a data file meant for 3D plot over a 2D grid. Inserts blank lines between lines of the grid.
# User must define {field} (1 or 2) with "-v". Assumes that the input file is sorthed by that field and
# then by the other field (2 or 1).

BEGIN { xold = ""; }

/[0-9]/ {
  x = $(field)
  if ((xold != "") && (x != xold)) { printf "\n"; }
  print $0;
  xold = x;
  next
}
