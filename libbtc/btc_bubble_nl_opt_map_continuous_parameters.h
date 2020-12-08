#ifndef btc_bubble_nl_opt_map_continuous_parameters_H
#define btc_bubble_nl_opt_map_continuous_parameters_H

/* Collecting adjustable integer paramters of a BTC price bubble model. */
/* Last edited on 2015-04-20 14:25:39 by stolfilocal */

void btc_bubble_nl_opt_map_continuous_parameters
  ( int npf, 
    double pf_lo[], 
    double pf[], 
    double pf_hi[],
    double x[]
  );
  /* Applies to each parameter {pf[ip]}, for {ip} in {0..npf-1}, an
    affine map that takes {[pf_lo[ip] _ pf_hi[ip]]} to {[-1 _ +1]}, and
    stores the result in {x[ip]}. The range must have non-zero width. */
    

#endif
