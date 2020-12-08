#! /usr/bin/gawk -f
# Last edited on 2017-06-15 18:06:49 by stolfilocal

BEGIN {
  # Prints the absolute and relative error in the conversion of JPEG's {JCS_RGB565}
  # color space samples, from {0..31} and {0..63} to {0..255} range.

    check_errors(5,31,8,255);
    check_errors(6,63,8,255);
    
    check_errors(5,31,16,65535);
    check_errors(6,63,16,65535);
    
}

function check_errors(ib,imax,ob,omax,   b,i,x,v,e,s,r,emax,smax,rmax)
{ 
  # Checks errors in sample conversion from {ib} bits ({0..imax}) 
  # to {ob} bits ({0..omax}).
  
  printf "----------------------------------------------------------------------\n" > "/dev/stderr";
  printf "conversion from %d bits (0..%d) to %d bits (0..%d)\n", ib,imax,ob,omax > "/dev/stderr";
  b = int(imax/2); # Rounding bias.
  printf "rounding bias = %d\n", b > "/dev/stderr";
  printf "\n" > "/dev/stderr";
  for (i = 1; i <= imax; i++)
    { x = i*omax/imax;              # "Exact" fractional converted value.
      v = int((i*omax + b)/imax); # Integer converted sample value.
      e = v - x;                    # Absolute error in sample value.
      s = e/omax;                   # Error relative to {omax}.
      r = e/x;                      # Error relative to exact value. 
      printf "%3d %12.6f %6d %+10.7f %+10.7f %+10.7f\n", i, x, v, e, s, r > "/dev/stderr";
      if (fabs(e) > emax) { emax = e; }
      if (fabs(s) > smax) { smax = s; }
      if (fabs(r) > rmax) { rmax = r; }
    }
    
  printf "\n" > "/dev/stderr";
  printf "max errors =  %+10.7f %+10.7f %+10.7f\n", emax, smax, rmax > "/dev/stderr";
  printf "----------------------------------------------------------------------\n" > "/dev/stderr";
}

function fabs(x) { if (x < 0) { return -x; } else { return x; } }
