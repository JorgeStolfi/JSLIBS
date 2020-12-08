#ifndef btc_bubble_fit_lsq_H
#define btc_bubble_fit_lsq_H

/* Least-squares fitting of the linear coefficients of a bubble model. */
/* Last edited on 2015-04-21 01:38:52 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_fit_lsq
  ( int nd, 
    char* dt[], 
    double ap[],
    double wt[],
    int nb, 
    btc_bubble_t bp[], 
    double bval[], 
    int maxIters,
    char* outPrefix
  );
  /* Adjusts the coefficients of the bubble functions for prices
    {ap[0..nd-1]} using robust least squares with {maxIters} iterations.
    The goal function is the RMS error between the modeled and the given
    price series, in linear scale, with each data poiny {ap[id]}
    weighted by {wt[id]}.
    
    The procedure assumes that the bubble function values are stored in
    {bval[0..nd*nb-1]}, with one row for each date. The fitted
    coefficients will be stored in in {bp[0..nb-1].coef}, ignoring their
    input values. The dates {dt[0..nd-1]} are used for debugging
    only. */

#endif
