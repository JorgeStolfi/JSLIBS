#ifndef btc_bubble_eval_rms_log_error_H
#define btc_bubble_eval_rms_log_error_H

/* Evaluating the rms log error of a BTC price bubble model. */
/* Last edited on 2015-04-22 20:44:53 by stolfilocal */

#include <btc_bubble_t.h>

double btc_bubble_eval_rms_log_error
  ( int nd, 
    double ap[],
    int id_ini,
    int id_fin,
    int nb,
    btc_bubble_t bp[],
    double bval[]
  ); 
  /* Returns the RMS difference between the decimal log of the observed
    price {ap[0..nd-1]} and the price defined by the BTC bubble model
    {bp[0..nb-1]}. Ignores dates when the price {ap[id]} is zero.
    Conisder only data points with indices in {id_in..id_fin}.
    Uses only the coefficients {.coef} from the {bp}
    parameters; assumes that the bubbles described by {bp} are evaluated
    in {bval}, smothed as appropriate. */

#endif
