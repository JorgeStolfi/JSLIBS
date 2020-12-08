#ifndef btc_bubble_nl_opt_set_continuous_variable_parameters_H
#define btc_bubble_nl_opt_set_continuous_variable_parameters_H

/* Setting the adjustable {double} parameters of a BTC price bubble model. */
/* Last edited on 2015-04-20 14:40:50 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_nl_opt_set_continuous_variable_parameters
  ( int npf, 
    double pf[], 
    int nb, 
    btc_bubble_t bp_lo[],
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[]
  );
  /* Sets the adjustable {double}-valued fields in the bubble parameters
    {bp[0..nb-1]} to the successive values {pf[0..npf-1]}. The
    adjustable parameters are identified by their non-trivial ranges in
    {bp_lo,bp_hi}, as in
    {btc_bubble_nl_opt_gather_continuous_variable_parameters}, and their
    number and order is the same as defined in that procedure. */

#endif
