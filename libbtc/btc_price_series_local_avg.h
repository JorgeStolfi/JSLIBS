#ifndef btc_price_series_local_avg_H
#define btc_price_series_local_avg_H

/* Local average of a BTC price series. */
/* Last edited on 2015-04-20 00:41:07 by stolfilocal */
  
double btc_price_series_local_avg(int nd, double val[], int id, int hrad);
  /* Assumes that {val[0..nd-1]} is the series to be smoothed.
    Computes the (log-scale) average of {val[kd]}  in a soft window 
    around day with index {id}. The window spans {2*hrad+1} days.
    Assumes that {val[kd]} is missing if it is zero.
    If {hrad} is zero, returns {val[id]}. */


#endif
