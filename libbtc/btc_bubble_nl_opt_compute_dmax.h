#ifndef btc_bubble_nl_opt_compute_dmax_H
#define btc_bubble_nl_opt_compute_dmax_H

/* Collecting adjustable integer paramters of a BTC price bubble model. */
/* Last edited on 2015-04-20 17:15:40 by stolfilocal */

double btc_bubble_nl_opt_compute_dmax(int npf, double x_ini[]);
  /* Computes the maximum distance from the vector {x[0..npf-1]} to 
    any point of the cube {[-1 _ +1]^npf}. */
    

#endif
