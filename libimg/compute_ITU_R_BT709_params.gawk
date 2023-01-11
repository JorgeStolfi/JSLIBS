#! /usr/bin/gawk -f
# Last edited on 2007-11-10 22:48:49 by stolfi

BEGIN{
  # Computes the correct values of the ITU-R Recommendation BT.709 
  # luminance encoding function.
  
  gamma = 0.450;         # Exponent; assumed correct.
  bmin = 0.099296826809; # Min value of {b} to try.
  bmax = 0.099296826810; # Max value of {b} to try.
  nb = 10;               # Number of steps in the {b} parameter.
  for (i = 0; i <= nb; i++)
    { s = (i - 0.0)/(nb - 0.0);
      b = (1 - s)*bmin + s*bmax;
      print_params(gamma, b);
    }
}

function print_params(g,b,  s,r,t,a,dfr,f0,f1,v0,v1)
{ # Prints parameters for exponent {g} and vertical offset {b}.
  s = b/(1+b)/(1-g); # Value of {r^g}.
  r = exp(log(s)/g); # Break argument value.
  dfr = (1+b)*g*s/r; # Slope of power curve at {v=r}.
  a = dfr;           # Slope of linear segment.
  f0 = a*r;                      # Function value just before break. 
  f1 = (1+b)*exp(g*log(r)) - b;  # Function value just after break.
  t = f0+r;                      # Break arg/val sum.
  v0 = s/a;                      # Inverse function value just before break.
  v1 = exp(log((r+b)/(1+b))/g);  # Inverse function value just after break.
  printf "           %23s  %23s\n", "--- encoding params ---", "--- decoding params ---"
  printf "offset     b =  %18.12f  B =  %18.12f\n",  b, b;  
  printf "exponent   g =  %18.12f  G =  %18.12f\n",  g, 1/g;
  printf "break arg  r =  %18.12f  R =  %18.12f\n",  r, s;  
  printf "break fun  s =  %18.12f  S =  %18.12f\n",  s, r;  
  printf "break sum  t =  %18.12f  T =  %18.12f\n",  t, t;  
  printf "init slope a =  %18.12f  A =  %18.12f\n",  a, 1/a;
  printf "left lim   f0 = %18.12f  F0 = %18.12f\n", f0, v0; 
  printf "right lim  f1 = %18.12f  F1 = %18.12f\n", f1, v1; 
  printf "           %23s  %23s\n", "-----------------------", "-----------------------"
  printf "\n"
}
