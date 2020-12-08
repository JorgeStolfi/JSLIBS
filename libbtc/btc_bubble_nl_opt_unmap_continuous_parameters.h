#ifndef btc_bubble_nl_opt_unmap_continuous_parameters_H
#define btc_bubble_nl_opt_unmap_continuous_parameters_H

/* Collecting adjustable integer paramters of a BTC price bubble model. */
/* Last edited on 2015-04-20 14:32:39 by stolfilocal */

double btc_bubble_nl_opt_unmap_continuous_parameters
  ( int npf,
    double x[], 
    double pf_lo[], 
    double pf[], 
    double pf_hi[]
  );
  /* Applies to each element {x[ip]}, for {ip} in {0..npf-1} an affine map 
    takes {[-1 _ +1]} to {[pf_lo[ip] _ pf_hi[ip]]}, and
    stores the result in {pf[ip]}, clipping it to that range (that must 
    have non-zero width).  Returns zero if {x} is inside the 
    cube {[-1 _ +1]^npf}, otherwise returns the Euclidean distance squared
    from {x} to that cube. */
    

#endif
