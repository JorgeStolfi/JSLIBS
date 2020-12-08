#ifndef btc_bubble_parms_copy_H
#define btc_bubble_parms_copy_H

/* Copying BTC price bubble parameters. */
/* Last edited on 2015-04-20 00:31:39 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_parms_copy(int nb, btc_bubble_t bp_a[], btc_bubble_t bp_b[]);
  /* Copies the bubble parameters {bp_a[0..nb-1]} to {bp_b[0..nb-1]}. */

#endif
