#! /usr/bin/python3
# Last edited on 2024-01-04 11:27:22 by stolfi

from math import exp, log, sqrt, pi;
from sys import stderr as err;

def gam(f, sig):
  return exp(-0.5*(f/sig)**2);
  
def W(f, s0, s1):
  return (1 - gam(f,s0))*gam(f, s1);
  
def app(f, s0, s1):
  return f*f*(1/s0**2 - 0.5*s1**2) - 1.5
  
fmax = 100
sig_lo = 20;
sig_hi = 70;
nf = 200;
for i in range(nf):
  f = fmax*i/nf;
  fp = f + 0.01;  wp = log(W(fp, sig_lo, sig_hi));
  fm = f - 0.01;  wm = log(W(fm, sig_lo, sig_hi));
  dw_num = (wp - wm)/(fp - fm);
  dw_app = app(f, sig_lo, sig_hi);
  err.write("f = %6.1f  num = %16.12f app = %16.12f\n" % (f,dw_num, dw_app));
  
  
