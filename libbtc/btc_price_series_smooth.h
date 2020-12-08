#ifndef btc_price_series_smooth_H
#define btc_price_series_smooth_H

/* Smoothing a BTC price series. */
/* Last edited on 2015-04-20 00:44:23 by stolfilocal */

void btc_price_series_smooth(int nd, double vi[], int hrad, double vo[]);
  /* Computes {vo[id]} as the average of {vi} in a window around {id},
    for {id} in {0..nd-1}. Assumes that {vi[kd]} is missing if it is zero.
    Sets {vo[id]} to zero if there is not enough data.
    If {hrad} is zero, sets {vo[id] = vi[id]} for all {id}. */


#endif
