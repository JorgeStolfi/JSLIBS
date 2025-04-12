#! /usr/bin/python3
# Last edited on 2025-03-16 18:02:33 by stolfi

import sys, os
from sys import stderr as err

def main():
  nz = 10;
  nw = 3;
  for kzmax in range(nz+2):
    Zmax = (kzmax-0.5)/nz;
    for kw in range(nw):
      W = (kw+1)/(nw + 2);
      fname = f"zwtest/Z{kzmax:02d}-W{kw:02d}.txt"
      wr = open(fname, "w")
      for kz in range(nz+1):
        Z = kz/nz;
        ef = (Z - Zmax)/W;
        vf = 1 + ef*ef
        f = 1/vf
        wr.write(f"{Z:6.3f} {f:8.5f}\n")
      wr.close()
  err.write("done.\n")

main()
