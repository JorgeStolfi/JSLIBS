#ifndef btc_bubble_nl_opt_check_bubble_parms_in_range_H
#define btc_bubble_nl_opt_check_bubble_parms_in_range_H

/* Checking whether guessed parameters are in specified range. */
/* Last edited on 2015-04-20 01:24:29 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_nl_opt_check_bubble_parms_in_range
  ( int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[]
  );
  /* Chechs whether all adjustable bubble parameters in {bp[0..nb-1]}
    are in the range defined by the corresponding parameters
    of {{bp_lo,bp_hi}[0..nb-1]}. Also checks whether the non-adjustable
    parameters have trivial ranges equal to the value in {bp}. */

#endif
