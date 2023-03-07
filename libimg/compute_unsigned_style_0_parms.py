#! /usr/bin/python3
# Last edited on 2023-03-05 21:24:53 by stolfi

# Computes the parameter A for 
# {frgb_path_map_unsigned_0}
#
# The path has 
#   {R(z) = aR*T(z)*sin(fR*2*PI*z) + z}
#   {B(z) = aB*T(z)*sin(fB*2*PI*z) + z}
# where {T(z) = S(z)*(2-S(z))} and {S(z) = 0.5*(1-cos(2*pi*z))}.
# The challenge is find max {aR,aB},
# so that {0 <= Y(R(z),0,B(z)) <= z - Y(0,1,0)}

from math import sin,cos,ceil,floor,sqrt,inf,pi
import sys

def main():
  stages = 4
  compute_aR_aB_table(stages)
  plot_some_curve()
  
def plot_some_curve():
  fR = 3; aR = 0.40;
  fB = 5; aB = 0.30;
  plot(fR,aR, fB,aB)
  return 0
  
def compute_aR_aB_table(stages):
  Cmax = 20
  print(Cmax)
  with open("out/uv.txt", "w") as f:
    for kc in range(2*Cmax+1):
      cycles = kc - Cmax;
      fR,fB = choose_freqs_close(cycles)
      assert check(fR,0,fB,0,4000), "duh?"
      aR,aB = find_aR_aB_hard(cycles,fR,fB,stages)
      f.write("%+3d   %+4d  %12.8f  %+4d %12.8f\n" % (cycles, fR, aR, fB, aB))
    
def choose_freqs_close(cycles):
  # Chooses two odd frequencies {fR,fB} that are relatively prime 
  # and close to each other.
  # If {cycles} is positive, {fR} will be {2*cycles+1} and
  # {fB} will be {fR-2 = 2*cycles-1}.  Since {fR} is odd,
  # the two will be relatively prime.
  # If {cycles} is negative, reverses {fR} and {fB}.
  # If {cycles} is zero, return {0,0}.
  if cycles == 0:
    fR = 0; fB = 0; 
  else:
    fR = 2*abs(cycles)+1;
    fB = 2*abs(cycles)-1;
  if cycles > 0:
    return fR,fB
  else:
    return fB,fR
    
def choose_freqs_gold(cycles):
  # Chooses two odd frequencies {fR,fB} that are relatively
  # prime, roughly in golden ratio.
  # If {cycles} is positive, {fR} will be {2*cycles+1} and
  # {fB} will be close to {phi*fR} where {phi} is the golden ratio.
  # If {cycles} is negative, reverses {fR} and {fB}.
  # If {cycles} is zero, return {0,0}.
  if cycles == 0:
    fR = 0; fB = 0; 
  elif abs(cycles) == 1:
    fR = 3; fB = 1;
  else:
    fR = 2*abs(cycles)+1;
    phi = (sqrt(5)-1)/2
    fB = int(floor(phi*fR));
    fB = fB + 1 - (fB%2);
    while fB == fR or fB <= 0 or gcd(fR,fB) != 1:
      fB = fB + 2
  if cycles > 0:
    return fR,fB
  else:
    return fB,fR
  
def plot(fR,aR,fB,aB):
  with open("out/plot.txt", "w") as f:
    nz = 200
    for kz in range(nz+1):
      z = kz/nz
      Rz,Bz,YRBz = RBY(fR,aR,fB,aB,z)
      Gz,Yz = GY(Rz,Bz,YRBz,z)
      Rmax = fmin(1.0, 
      f.write("%4d %12.8f" % (kz, z));  
      f.write("  %12.8f %12.8f %12.8f" % (Rz,Gz,Bz));  
      f.write("  %12.8f %12.8f\n" % (Yz,YRBz));  
      f.write("  %12.8f %12.8f\n" % (Yz,YRBz));  
    
def find_aR_aB_hard(cycles,fR,fB,stages):
  na = 20;
  nz = 250;
  stages = 4;
  aR = 0.5; aB = 0.5; d = 0.5;
  for ks in range(stages):
    aR, aB = find_aR_aB_aux(fR, fB, aR-d,aR+d, aB-d,aB+d, na, nz);
    d = d/10
    nz = 4*nz
  return aR, aB

def find_aR_aB_fiat_gold(cycles,fR,fB,stages):
  if cycles == 0:
    assert fR == 0 and fB == 0
    aR = 0.000; aB = 0.000
  elif cycles == -1: 
    assert fR == 1 and fB == 3
    aR = 0.383; aB = 0.663
  elif cycles == -2:
    assert fR == 3 and fB == 5
    aR = 0.661; aB = 0.120
  elif cycles == -3:
    assert fR == 5 and fB == 7
    aR = 0.396; aB = 0.336
  elif cycles == -4:
    fR == 5 and fB == 9
    aR = 0.396; aB = 0.338
  elif cycles == +1:
    assert fR == 3 and fB == 1
    aR = 0.655; aB = 0.335
  elif cycles == +2:
    assert fR == 5 and fB == 3
    aR = 0.396; aB = 0.663
  elif cycles == +3:
    fR == 7 and fB == 5
    aR = 0.336; aB = 0.396
  elif cycles == +4:
    fR == 9 and fB == 5
    aR = 0.338; aB = 0.396
  else:
    aR, aB = find_aR_aB_hard(cycles,fR,fB,stages)
  nz = 10000;
  check(fR,aR,fB,aB,nz) 
  return aR, aB
  
def find_aR_aB_aux(fR,fB, aRlo,aRhi, aBlo,aBhi, na, nz):
  sys.stderr.write("%+4d %+4d  %12.8f %12.8f  %12.8f %12.8f  %6d %6d\n" % (fR,fB,aRlo,aRhi,aBlo,aBhi,na,nz))
  aRbest = -inf
  aBbest = -inf
  for kR in range(na+1):
    aR = aRlo + kR*(aRhi - aRlo)/na
    for kB in range(na+1):
      aB = aBlo + kB*(aBhi - aBlo)/na
      ok = check(fR,aR,fB,aB,nz)
      if ok and (aR > aRbest or (aR == aRbest and aB > aBbest)):
        aRbest = aR; aBbest = aB
  return aRbest, aBbest
  
def check(fR,aR,fB,aB,nz):
  # sys.stderr.write("%+4d %12.8f %+4d %12.8f %6d\n" % (fR,aR,fB,aB,nz))
  if aR < 0 or aB < 0:
    return False
  for kz in range(nz+1):
    z = kz/nz
    Rz,Bz,YRBz = RBY(fR,aR,fB,aB,z)
    if Rz < 0 or Rz > 1: return False
    if Bz < 0 or Bz > 1: return False
    Gz,Yz = GY(Rz,Bz,YRBz,z)
    if Gz < 0 or Gz > 1: return False
  return True
  
def GY(Rz,Bz,YRBz,z):
  Gz = (z - YRBz)/0.5866116
  Yz = 0.298911*Rz +0.5866116*Gz +0.114478*Bz
  return Gz,Yz

def RBY(fR,aR,fB,aB,z):
  Sz = 0.5*(1-cos(2*pi*z))
  Tz = Sz*(2 - Sz)
  Rz = z + aR*Tz*0.5*(1-cos(fR*pi*z))
  Bz = z + aB*Tz*0.5*(1-cos(fB*pi*z))
  YRBz = 0.298911*Rz +0.114478*Bz
  return Rz,Bz,YRBz

def gcd(x,y):
  x = abs(x)
  y = abs(y)
  while y != 0:
    r = x % y; x = y; y = r
  return x

main()
