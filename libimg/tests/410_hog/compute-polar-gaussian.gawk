#! /usr/bin/gawk -f
# Last edited on 2012-01-10 00:02:00 by stolfilocal

# Computes the approximate angular distribution of a two-dimensional
# Gaussian distribution and compares with an approximate formula.

BEGIN {
  if (na == "") { arg_error(("must define {na}")); }
  nd = 6;          # Steps in center distance {d}.
  pi = 3.1415926;  # The Archimedean constant.
  split("", hist); # Histogram, indexed by {id} and {ia}. 
  dmin = 0.25;
  d = dmin;
  for (id = 0; id < nd; id++) {
    printf "d%d = %10.6f\n", id, d > "/dev/stderr";
    compute_azimuth_hist(id,d,na,hist);
    d = 2*d;
  }

  for (ia = 0; ia <= na; ia++) {
    fa = ia/na; # Bin central angle
    a = pi*fa;
    printf "%14.10f %14.10f ", fa, a;
    for (id = 0; id < nd; id++) {
      printf " %24.16e", hist[id,ia];
    }
    printf "\n";
    printf "angle = %10.6f half-turns = %10.6f radians\n", fa, a > "/dev/stderr";
  }
}

function compute_azimuth_hist(id,d,na,hist,   ia,nj,ja,fj,a,sum,binsum   )
  {
    # Assumes that {(x,y)} has a two-dimensional Gaussian distribution {Pr(x,y)} with
    # mean {(d,0)} and standard deviation 1. Stores in {hist[id,0..na]} the distribution of the 
    # azimuth of {(x,y)}, sampled at {na} angles between 0 and {PI}, inclusive.
    # The distribution is normalized to unit sum.
    
    sum = 0;
    nj = 1;
    for (ia = 0; ia <= na; ia++) {
      # Sample the distribution over {[pi*(ia-1)/na _ pi*(ia+1)/na]} with Hann weight.
      binsum = 0;
      wsum = 0;
      for (ja = 1-nj; ja <= nj-1; ja++) {
        fj = ja/nj;
        a = pi*(ia + fj)/na;
        w = 0.5*(1 + cos(pi*fj));  # Hann bin weight.
        binsum += w*radial_prob(d,a);
        wsum += w;
      }
      binsum = binsum/wsum;
      hist[id,ia] = binsum;
      sum += binsum;
    }

    # Normalize to unit integral:
    for (ia = 0; ia <= na; ia++) {
      hist[id,ia] /= (sum/na);
    }
  }

function radial_prob(d,a,   ca,sa,b,emax,r0,hr,rmin,rmax,nr,ir,r,x,y,pr,sum)
  {
    # Assumes that {(x,y)} has a two-dimensional Gaussian distribution {Pr(x,y)} with
    # mean {(d,0)} and standard deviation 1. Computes the (un-normalized) probability
    # density of the azimuth of {(x,y)} being approximately {a}.

    ca = cos(a);
    sa = sin(a); 
    # Azimuth distribution is symmetric around zero:
    if (sa < 0) { sa = -sa; }

    emax = 10;   # Assume {Pr(x,y)=0} if {dist((x,y),(d,0)) > emax}.
    b = sa*d;   # Distance from center to ray.
    if (emax <= b) {
      # Ray does not intercept {Pr}:
      return 0;
    } else {
      r0 = ca*d;                       # Radius {r} closest to center.
      hr = sqrt(emax*emax - b*b);      # Half of ray segment that is inside {Pr}.
      rmin = (hr >= r0 ? 0 : r0 - hr); # Start of integration interval.
      rmax = r0 + hr;                  # End of integration interval.
      nr = 600;   # Radial integration steps.
      sum = 0;
      for (ir = 0; ir < nr; ir++) {
        r = ((nr - ir - 0.5)*rmin + (ir + 0.5)*rmax)/nr;
        x = r*ca;
        y = r*sa;
        pr = exp(-0.5*((x-d)*(x-d) + y*y));
        sum = sum + pr*r;
      }
      return (rmax - rmin)*sum/nr;
    }
  }
