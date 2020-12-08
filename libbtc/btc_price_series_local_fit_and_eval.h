#ifndef btc_price_series_local_fit_and_eval_H
#define btc_price_series_local_fit_and_eval_H

/* Local polynomial smoothing/interpolation of a BTC price series. */
/* Last edited on 2015-04-20 00:39:27 by stolfilocal */
    
double btc_price_series_local_fit_and_eval(int deg, int hrad, double val[], double wht[]);
  /* Assumes that {val[0..nw-1]} and {wht[0..nw-1]} are 
    the data and weights within some window of size {nw = 2*hrad+1}
    Fits a polynomial {P(j)} of degree {deg} to the data points {(j,log(val[hrad+j]))} with
    weight {wht[hrad+j]}, for {j} in {-hrad..+hrad}.  Then evaluates {exp(P(0))}.
    Returns 0.0 if there is not enough data to fit the polynomial. */


#endif
