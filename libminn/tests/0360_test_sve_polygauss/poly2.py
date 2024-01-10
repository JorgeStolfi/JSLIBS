#! /usr/bin/python3
# Last edited on 2024-01-08 05:25:46 by stolfi

from math import exp, log, sqrt;
from sys import stdout as out, stderr as err;

def polygauss(x,np,sigma):
  F = 0
  for ip in range(np):
    avg = ip/(np-1)
    F += gauss(x,avg,sigma)
  return F
  
def gauss(x,avg,sigma):
  z = (x - avg)/sigma
  return exp(-0.5*z*z)
  
def plot():
  np = 3
  nx = 200
  hs = 3
  sigma0 = 0.5/(np-1) + 0.05*(np-2);
  for ix in range(nx+1):
    x = ix/nx
    out.write(f"{x:12.9f} ");
    for ks in range(2*hs+1):
      sigma = sigma0 + 0.05*(ks - hs)/hs/(np-1);
      if ix == 0:
        err.write(f"  {x:12.9f} {sigma:7.4f}\n")
      f = polygauss(x,np,sigma);
      out.write(f" {f:16.12f}");
    out.write("\n")
  out.flush()
  
plot()
