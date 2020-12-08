/* See {btc_bubble_parms_copy.h} */
/* Last edited on 2015-04-20 00:30:00 by stolfilocal */

#define _GNU_SOURCE

#include <btc_bubble_t.h>
#include <btc_bubble_parms_copy.h>

void btc_bubble_parms_copy(int nb, btc_bubble_t bp_a[], btc_bubble_t bp_b[])
  { 
    int ib;
    for (ib = 0; ib < nb; ib++) { bp_b[ib] = bp_a[ib]; }
  }

